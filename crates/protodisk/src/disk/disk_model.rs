//! Trait abstraction for protoplanetary disk models.
//!
//! This trait enables both analytical power-law disks (for fast statistical
//! generation) and grid-based evolving disks (for emergent simulation) to
//! share physics calculations.
//!
//! # Design
//!
//! The trait separates "what defines a disk" from "what we can compute from it":
//!
//! | Category | Methods | Rationale |
//! |----------|---------|-----------|
//! | **Required** | 6 methods | Minimal set that varies between representations |
//! | **Override-able** | `pressure_gradient_log` | Numerical default; analytical override for power-law |
//! | **Derived** | 15+ methods | Physics identical regardless of representation |
//!
//! A new disk type only needs to implement 6 methods to get full functionality.

use units::{
    AngularVelocity, Density, Length, Mass, Pressure, SurfaceDensity, Temperature, Time, Velocity,
};

use crate::disk::constants::{G, K_B, M_PROTON, MU, PI};

/// Molecular cross section for H2 (cm²)
const SIGMA_MOL: f64 = 2e-15;

/// A protoplanetary disk model.
///
/// Implementors provide the fundamental disk properties (surface density,
/// temperature, stellar mass, etc.) and get derived quantities (sound speed,
/// scale height, pressure, viscosity, etc.) for free via default implementations.
pub trait DiskModel {
    // =========================================================================
    // Required methods - these define the disk
    // =========================================================================

    /// Surface density at radius r.
    fn surface_density(&self, r: Length) -> SurfaceDensity;

    /// Temperature at radius r.
    fn temperature(&self, r: Length) -> Temperature;

    /// Central stellar mass.
    fn stellar_mass(&self) -> Mass;

    /// Shakura-Sunyaev viscosity parameter α.
    fn alpha(&self) -> f64;

    /// Inner edge of the disk.
    fn inner_radius(&self) -> Length;

    /// Outer edge of the disk.
    fn outer_radius(&self) -> Length;

    // =========================================================================
    // Override-able methods - have sensible defaults
    // =========================================================================

    /// Logarithmic pressure gradient: d ln P / d ln r
    ///
    /// Default implementation uses numerical central difference.
    /// Power-law disks should override with the analytical expression.
    fn pressure_gradient_log(&self, r: Length) -> f64 {
        // Use central difference in log space
        let delta = 0.01; // 1% perturbation
        let r_cm = r.to_cm();

        let r_minus = Length::from_cm(r_cm * (1.0 - delta));
        let r_plus = Length::from_cm(r_cm * (1.0 + delta));

        let p_minus = self.pressure(r_minus).to_dyn_per_cm2();
        let p_plus = self.pressure(r_plus).to_dyn_per_cm2();

        // d ln P / d ln r ≈ (ln P+ - ln P-) / (ln r+ - ln r-)
        let ln_p_diff = p_plus.ln() - p_minus.ln();
        let ln_r_diff = (r_cm * (1.0 + delta)).ln() - (r_cm * (1.0 - delta)).ln();

        ln_p_diff / ln_r_diff
    }

    // =========================================================================
    // Derived methods - computed from the required methods
    // =========================================================================

    /// Keplerian orbital frequency at radius r.
    /// Ω_K = √(GM_*/r³)
    fn orbital_frequency(&self, r: Length) -> AngularVelocity {
        let r_cm = r.to_cm();
        let omega = (G * self.stellar_mass().to_grams() / r_cm.powi(3)).sqrt();
        AngularVelocity::from_rad_per_sec(omega)
    }

    /// Keplerian orbital velocity at radius r.
    /// v_K = √(GM_*/r)
    fn keplerian_velocity(&self, r: Length) -> Velocity {
        let r_cm = r.to_cm();
        let v_k = (G * self.stellar_mass().to_grams() / r_cm).sqrt();
        Velocity::from_cm_per_sec(v_k)
    }

    /// Orbital period at radius r.
    /// P = 2π/Ω_K
    fn orbital_period(&self, r: Length) -> Time {
        let omega = self.orbital_frequency(r);
        Time::from_seconds(2.0 * PI / omega.to_rad_per_sec())
    }

    /// Isothermal sound speed at radius r.
    /// c_s = √(k_B T / (μ m_p))
    fn sound_speed(&self, r: Length) -> Velocity {
        let t = self.temperature(r).to_kelvin();
        let c_s = (K_B * t / (MU * M_PROTON)).sqrt();
        Velocity::from_cm_per_sec(c_s)
    }

    /// Thermal velocity (mean molecular speed) at radius r.
    /// v_th = √(8/π) × c_s ≈ 1.60 × c_s
    fn thermal_velocity(&self, r: Length) -> Velocity {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        Velocity::from_cm_per_sec(c_s * (8.0 / PI).sqrt())
    }

    /// Disk scale height at radius r.
    /// h = c_s / Ω_K
    fn scale_height(&self, r: Length) -> Length {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let omega = self.orbital_frequency(r).to_rad_per_sec();
        Length::from_cm(c_s / omega)
    }

    /// Disk aspect ratio at radius r.
    /// h/r = c_s / v_K
    fn aspect_ratio(&self, r: Length) -> f64 {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let v_k = self.keplerian_velocity(r).to_cm_per_sec();
        c_s / v_k
    }

    /// Midplane gas density at radius r.
    /// ρ = Σ / (√(2π) × h)
    ///
    /// Assumes a Gaussian vertical density profile.
    fn midplane_density(&self, r: Length) -> Density {
        let sigma = self.surface_density(r).to_grams_per_cm2();
        let h = self.scale_height(r).to_cm();
        let rho = sigma / ((2.0 * PI).sqrt() * h);
        Density::from_grams_per_cm3(rho)
    }

    /// Midplane pressure at radius r.
    /// P = ρ × c_s²
    fn pressure(&self, r: Length) -> Pressure {
        let rho = self.midplane_density(r).to_grams_per_cm3();
        let c_s = self.sound_speed(r).to_cm_per_sec();
        Pressure::from_dyn_per_cm2(rho * c_s.powi(2))
    }

    /// Pressure gradient parameter η.
    ///
    /// η = -(h/r)² × (1/2) × d ln P / d ln r
    ///
    /// This controls sub-Keplerian gas rotation: v_φ,gas = v_K × (1 - η)
    /// Typical values: η ≈ 0.002-0.005
    fn pressure_gradient_parameter(&self, r: Length) -> f64 {
        let h_over_r = self.aspect_ratio(r);
        let d_ln_p = self.pressure_gradient_log(r);
        -h_over_r.powi(2) * 0.5 * d_ln_p
    }

    /// Gas velocity relative to Keplerian.
    /// Δv = v_K - v_φ,gas = η × v_K
    ///
    /// This is the "headwind" experienced by solid particles.
    fn sub_keplerian_velocity(&self, r: Length) -> Velocity {
        let eta = self.pressure_gradient_parameter(r);
        let v_k = self.keplerian_velocity(r).to_cm_per_sec();
        Velocity::from_cm_per_sec(eta * v_k)
    }

    /// Kinematic viscosity at radius r using α-prescription.
    /// ν = α × c_s × h
    fn viscosity(&self, r: Length) -> f64 {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let h = self.scale_height(r).to_cm();
        self.alpha() * c_s * h
    }

    /// Viscous timescale at radius r.
    /// t_visc = r² / ν
    fn viscous_timescale(&self, r: Length) -> Time {
        let r_cm = r.to_cm();
        let nu = self.viscosity(r);
        Time::from_seconds(r_cm.powi(2) / nu)
    }

    /// Mean free path of gas molecules at radius r.
    /// λ = μ m_p / (σ_mol × ρ)
    fn mean_free_path(&self, r: Length) -> Length {
        let rho = self.midplane_density(r).to_grams_per_cm3();
        let lambda = MU * M_PROTON / (SIGMA_MOL * rho);
        Length::from_cm(lambda)
    }

    /// Check if radius is within disk bounds.
    fn is_valid_radius(&self, r: Length) -> bool {
        let r_cm = r.to_cm();
        r_cm >= self.inner_radius().to_cm() && r_cm <= self.outer_radius().to_cm()
    }
}

/// Trait for disk models that can compute total mass.
///
/// This is separate from `DiskModel` because the integration method
/// differs between representations:
/// - Power-law disks use closed-form integrals
/// - Grid disks use numerical integration
pub trait DiskMass: DiskModel {
    /// Total disk mass between inner and outer radius.
    fn total_mass(&self) -> Mass;
}
