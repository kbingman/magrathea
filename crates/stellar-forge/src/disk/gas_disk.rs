//! Gas disk structure and evolution
//!
//! Provides a protoplanetary gas disk with power-law surface density and
//! temperature profiles. All other quantities (scale height, midplane density,
//! pressure, etc.) are derived from these.
//!
//! # Physics
//!
//! The disk is assumed to be vertically isothermal and in hydrostatic equilibrium.
//! Key relations:
//!
//! - Scale height: h = c_s / Ω_K
//! - Midplane density: ρ = Σ / (√(2π) h)
//! - Pressure: P = ρ c_s²
//! - Pressure gradient parameter: η = -(h/r)² × (1/2) × d ln P / d ln r
//!
//! The pressure gradient parameter η controls the sub-Keplerian rotation of
//! gas and drives radial drift of solid particles.

use serde::{Deserialize, Serialize};
use units::{
    AngularVelocity, Density, Length, Mass, Pressure, SurfaceDensity, Temperature, Time, Velocity,
};

use crate::{
    MainSequenceStar,
    disk::constants::{G, K_B, M_PROTON, MU, PI},
    solar_analog,
};

/// A protoplanetary gas disk with power-law profiles.
///
/// Surface density: Σ(r) = Σ_0 × (r / r_0)^(-p)
/// Temperature: T(r) = T_0 × (r / r_0)^(-q)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GasDisk {
    /// Inner edge of the disk
    pub inner_radius: Length,
    /// Outer edge of the disk
    pub outer_radius: Length,

    /// Reference radius (typically 1 AU)
    pub r_0: Length,
    /// Surface density at r_0
    pub sigma_0: SurfaceDensity,
    /// Surface density power-law exponent (Σ ∝ r^(-p))
    /// Typical value: 1.0 (MMSN) to 1.5
    pub sigma_exponent: f64,

    /// Temperature at r_0
    pub temperature_0: Temperature,
    /// Temperature power-law exponent (T ∝ r^(-q))
    /// Typical value: 0.5 (flared disk) to 0.75 (flat disk)
    pub temp_exponent: f64,

    /// Central stellar mass
    pub stellar_mass: Mass,

    /// Shakura-Sunyaev viscosity parameter
    /// Typical value: 1e-3 to 1e-2
    pub alpha: f64,
}

impl GasDisk {
    /// Create a disk scaled to stellar properties.
    ///
    /// Temperature scaling for a passively irradiated disk:
    /// T(r) ∝ L*^(1/4) × r^(-1/2)
    ///
    /// At 1 AU around the Sun, T ≈ 280K. For other stars:
    /// T_0 = 280K × (L/L☉)^(1/4)
    ///
    /// Surface density is scaled to maintain similar disk-to-star mass ratio:
    /// Σ_0 = 1700 g/cm² × (M/M☉)
    pub fn for_star(star: &MainSequenceStar) -> Self {
        let mass_solar = star.mass.to_solar_masses();
        let luminosity_solar = star.luminosity; // already in solar luminosities

        // Temperature at 1 AU scales as L^(1/4)
        // Solar value: 280K at 1 AU
        let temp_at_1au = 280.0 * luminosity_solar.powf(0.25);

        // Surface density scales with stellar mass to maintain
        // similar disk-to-star mass ratio
        // Solar value: ~1700 g/cm² at 1 AU
        let sigma_at_1au = 1700.0 * mass_solar;

        // Disk extent scales with stellar mass (approximation)
        // More massive stars have more extended disks
        let inner = 0.1 * mass_solar.powf(0.5);
        let outer = 100.0 * mass_solar.powf(0.5);

        Self {
            inner_radius: Length::from_au(inner),
            outer_radius: Length::from_au(outer),
            r_0: Length::from_au(1.0),
            sigma_0: SurfaceDensity::from_grams_per_cm2(sigma_at_1au),
            sigma_exponent: 1.0,
            temperature_0: Temperature::from_kelvin(temp_at_1au),
            temp_exponent: 0.5,
            stellar_mass: star.mass,
            alpha: 1e-3,
        }
    }

    /// Create a Minimum Mass Solar Nebula disk around a solar-mass star.
    ///
    /// MMSN parameters from Hayashi (1981):
    /// - Σ = 1700 × (r/AU)^(-1.5) g/cm² (gas)
    /// - T = 280 × (r/AU)^(-0.5) K
    ///
    /// Note: We use a slightly shallower Σ exponent (1.0) which is more
    /// consistent with observed disks.
    pub fn mmsn() -> Self {
        Self::for_star(&solar_analog())
    }

    /// Create a disk with custom parameters.
    pub fn new(
        stellar_mass: Mass,
        sigma_0: SurfaceDensity,
        sigma_exponent: f64,
        temperature_0: Temperature,
        temp_exponent: f64,
        alpha: f64,
    ) -> Self {
        Self {
            inner_radius: Length::from_au(0.1),
            outer_radius: Length::from_au(100.0),
            r_0: Length::from_au(1.0),
            sigma_0,
            sigma_exponent,
            temperature_0,
            temp_exponent,
            stellar_mass,
            alpha,
        }
    }

    // =========================================================================
    // Primary profiles (power laws)
    // =========================================================================

    /// Surface density at radius r.
    /// Σ(r) = Σ_0 × (r / r_0)^(-p)
    pub fn surface_density(&self, r: Length) -> SurfaceDensity {
        let ratio = r.to_cm() / self.r_0.to_cm();
        SurfaceDensity::from_grams_per_cm2(
            self.sigma_0.to_grams_per_cm2() * ratio.powf(-self.sigma_exponent),
        )
    }

    /// Temperature at radius r.
    /// T(r) = T_0 × (r / r_0)^(-q)
    pub fn temperature(&self, r: Length) -> Temperature {
        let ratio = r.to_cm() / self.r_0.to_cm();
        Temperature::from_kelvin(self.temperature_0.to_kelvin() * ratio.powf(-self.temp_exponent))
    }

    // =========================================================================
    // Derived quantities
    // =========================================================================

    /// Keplerian orbital frequency at radius r.
    /// Ω_K = √(G M_* / r³)
    pub fn orbital_frequency(&self, r: Length) -> AngularVelocity {
        let r_cm = r.to_cm();
        let omega = (G * self.stellar_mass.to_grams() / r_cm.powi(3)).sqrt();

        AngularVelocity::from_rad_per_sec(omega)
    }

    /// Keplerian orbital velocity at radius r.
    /// v_K = r × Ω_K = √(G M_* / r)
    pub fn keplerian_velocity(&self, r: Length) -> Velocity {
        let r_cm = r.to_cm();
        let v_k = (G * self.stellar_mass.to_grams() / r_cm).sqrt();

        Velocity::from_cm_per_sec(v_k)
    }

    /// Orbital period at radius r.
    pub fn orbital_period(&self, r: Length) -> Time {
        let omega = self.orbital_frequency(r);
        Time::from_seconds(2.0 * PI / omega.to_rad_per_sec())
    }

    /// Isothermal sound speed at radius r.
    /// c_s = √(k_B T / (μ m_p))
    pub fn sound_speed(&self, r: Length) -> Velocity {
        let t = self.temperature(r).to_kelvin();
        let c_s = (K_B * t / (MU * M_PROTON)).sqrt();
        Velocity::from_cm_per_sec(c_s)
    }

    /// Thermal velocity (mean molecular speed) at radius r.
    /// v_th = √(8/π) × c_s ≈ 1.60 × c_s
    pub fn thermal_velocity(&self, r: Length) -> Velocity {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        Velocity::from_cm_per_sec(c_s * (8.0 / PI).sqrt())
    }

    /// Disk scale height at radius r.
    /// h = c_s / Ω_K
    pub fn scale_height(&self, r: Length) -> Length {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let omega = self.orbital_frequency(r).to_rad_per_sec();

        Length::from_cm(c_s / omega)
    }

    /// Disk aspect ratio at radius r.
    /// h/r = c_s / v_K
    pub fn aspect_ratio(&self, r: Length) -> f64 {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let v_k = self.keplerian_velocity(r).to_cm_per_sec();

        c_s / v_k
    }

    /// Midplane gas density at radius r.
    /// ρ = Σ / (√(2π) × h)
    ///
    /// This assumes a Gaussian vertical density profile.
    pub fn midplane_density(&self, r: Length) -> Density {
        let sigma = self.surface_density(r).to_grams_per_cm2();
        let h = self.scale_height(r).to_cm();
        let rho = sigma / ((2.0 * PI).sqrt() * h);

        Density::from_grams_per_cm3(rho)
    }

    /// Midplane pressure at radius r.
    /// P = ρ × c_s²
    pub fn pressure(&self, r: Length) -> Pressure {
        let rho = self.midplane_density(r).to_grams_per_cm3();
        let c_s = self.sound_speed(r).to_cm_per_sec();

        Pressure::from_dyn_per_cm2(rho * c_s.powi(2))
    }

    /// Logarithmic pressure gradient: d ln P / d ln r
    ///
    /// For power-law profiles with Σ ∝ r^(-p) and T ∝ r^(-q):
    /// - ρ ∝ Σ/h ∝ r^(-p) / r^((3-q)/2) = r^(-(p + (3-q)/2))
    /// - c_s² ∝ T ∝ r^(-q)
    /// - P ∝ ρ c_s² ∝ r^(-(p + (3-q)/2 + q)) = r^(-(p + (3+q)/2))
    ///
    /// Therefore: d ln P / d ln r = -(p + (3+q)/2)
    pub fn pressure_gradient_log(&self, _r: Length) -> f64 {
        // For power-law disk, this is constant
        -(self.sigma_exponent + (3.0 + self.temp_exponent) / 2.0)
    }

    /// Pressure gradient parameter η.
    ///
    /// η = -(h/r)² × (1/2) × d ln P / d ln r
    ///
    /// This parameter controls the sub-Keplerian rotation of gas:
    /// v_φ,gas = v_K × (1 - η)
    ///
    /// For typical disk parameters, η ≈ 0.002-0.005.
    pub fn pressure_gradient_parameter(&self, r: Length) -> f64 {
        let h_over_r = self.aspect_ratio(r);
        let d_ln_p = self.pressure_gradient_log(r);
        -h_over_r.powi(2) * 0.5 * d_ln_p
    }

    /// Gas velocity relative to Keplerian.
    /// Δv = v_K - v_φ,gas = η × v_K
    ///
    /// This is the "headwind" experienced by solid particles on
    /// Keplerian orbits.
    pub fn sub_keplerian_velocity(&self, r: Length) -> Velocity {
        let eta = self.pressure_gradient_parameter(r);
        let v_k = self.keplerian_velocity(r).to_cm_per_sec();
        Velocity::from_cm_per_sec(eta * v_k)
    }

    /// Viscosity at radius r using α-prescription.
    /// ν = α × c_s × h
    pub fn viscosity(&self, r: Length) -> f64 {
        let c_s = self.sound_speed(r).to_cm_per_sec();
        let h = self.scale_height(r).to_cm();

        self.alpha * c_s * h
    }

    /// Viscous timescale at radius r.
    /// t_visc = r² / ν
    pub fn viscous_timescale(&self, r: Length) -> Time {
        let r_cm = r.to_cm();
        let nu = self.viscosity(r);
        Time::from_seconds(r_cm.powi(2) / nu)
    }

    /// Mean free path of gas molecules at radius r.
    /// λ = μ m_p / (σ_mol × ρ)
    ///
    /// Uses σ_mol ≈ 2×10^(-15) cm² for H2.
    pub fn mean_free_path(&self, r: Length) -> Length {
        const SIGMA_MOL: f64 = 2e-15; // cm²
        let rho = self.midplane_density(r).to_grams_per_cm3();
        let lambda = MU * M_PROTON / (SIGMA_MOL * rho);

        Length::from_cm(lambda)
    }

    /// Total disk mass between inner and outer radius.
    /// M_disk = ∫ 2πr Σ(r) dr
    ///
    /// For Σ ∝ r^(-p):
    /// - p ≠ 2: M = 2π Σ_0 r_0^p × (r_out^(2-p) - r_in^(2-p)) / (2-p)
    /// - p = 2: M = 2π Σ_0 r_0² × ln(r_out/r_in)
    pub fn total_mass(&self) -> Mass {
        let sigma_0 = self.sigma_0.to_grams_per_cm2();
        let r_0 = self.r_0.to_cm();
        let r_in = self.inner_radius.to_cm();
        let r_out = self.outer_radius.to_cm();
        let p = self.sigma_exponent;

        let mass = if (p - 2.0).abs() < 1e-10 {
            // p ≈ 2: logarithmic integral
            2.0 * PI * sigma_0 * r_0.powi(2) * (r_out / r_in).ln()
        } else {
            // General case
            let factor = 2.0 * PI * sigma_0 * r_0.powf(p) / (2.0 - p);
            factor * (r_out.powf(2.0 - p) - r_in.powf(2.0 - p))
        };

        Mass::from_grams(mass)
    }

    // =========================================================================
    // Validation helpers
    // =========================================================================

    /// Check if radius is within disk bounds.
    pub fn is_valid_radius(&self, r: Length) -> bool {
        r.to_cm() >= self.inner_radius.to_cm() && r.to_cm() <= self.outer_radius.to_cm()
    }
}
