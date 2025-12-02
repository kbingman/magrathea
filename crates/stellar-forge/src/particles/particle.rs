//! Particle properties and aerodynamic drag.
//!
//! Implements drag physics following Weidenschilling (1977) and
//! Birnstiel et al. (2010).

use units::{Density, Length, Time, Velocity};

use crate::disk::DiskModel;
use crate::disk::constants::PI;

/// Aerodynamic drag regime for a particle in gas.
///
/// The regime depends on the particle size relative to the mean free path
/// and the Reynolds number of the flow around the particle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DragRegime {
    /// Particle smaller than gas mean free path (s < 9λ/4).
    /// Most particles in protoplanetary disks are in this regime.
    /// Drag force: F_D = (4/3) π s² ρ_g v_th Δv
    Epstein,

    /// Particle larger than mean free path, laminar flow (Re < 1).
    /// Drag force: F_D = 6π μ s Δv (Stokes law)
    Stokes,

    /// Particle larger than mean free path, turbulent flow (Re > 1).
    /// More complex drag law; not yet implemented.
    Turbulent,
}

/// A solid particle in a protoplanetary disk.
///
/// Represents dust grains, pebbles, or boulders characterized by
/// their size and internal density.
#[derive(Debug, Clone)]
pub struct Particle {
    /// Particle radius (cm)
    size: f64,

    /// Internal material density (g/cm³)
    /// Typical values: ~1 g/cm³ (ice), ~3 g/cm³ (silicate)
    material_density: f64,
}

impl Particle {
    /// Create a new particle with given size and material density.
    ///
    /// # Arguments
    /// * `size` - Particle radius
    /// * `material_density` - Internal density of the particle material
    ///
    /// # Example
    /// ```
    /// use stellar_forge::particles::Particle;
    /// use units::{Length, Density};
    ///
    /// // 1 mm silicate grain
    /// let grain = Particle::new(
    ///     Length::from_cm(0.1),
    ///     Density::from_grams_per_cm3(3.0),
    /// );
    /// ```
    pub fn new(size: Length, material_density: Density) -> Self {
        Self {
            size: size.to_cm(),
            material_density: material_density.to_grams_per_cm3(),
        }
    }

    /// Particle radius.
    pub fn size(&self) -> Length {
        Length::from_cm(self.size)
    }

    /// Internal material density.
    pub fn material_density(&self) -> Density {
        Density::from_grams_per_cm3(self.material_density)
    }

    /// Particle mass assuming a sphere.
    /// m = (4/3) π s³ ρ_m
    pub fn mass(&self) -> f64 {
        (4.0 / 3.0) * PI * self.size.powi(3) * self.material_density
    }

    /// Determine the drag regime for this particle at radius r in the disk.
    ///
    /// The regime depends on particle size relative to the gas mean free path.
    /// For most disk particles (μm to cm), Epstein drag applies.
    pub fn drag_regime<D: DiskModel>(&self, disk: &D, r: Length) -> DragRegime {
        let mfp = disk.mean_free_path(r).to_cm();

        // Epstein-Stokes transition at s = 9λ/4
        if self.size < 9.0 * mfp / 4.0 {
            DragRegime::Epstein
        } else {
            // For now, assume Stokes regime for larger particles
            // A full implementation would check Reynolds number
            DragRegime::Stokes
        }
    }

    /// Stopping time (friction time) of the particle.
    ///
    /// The stopping time measures how quickly the particle velocity
    /// relaxes to the gas velocity due to drag.
    ///
    /// # Epstein regime (s < 9λ/4)
    /// t_s = (ρ_m × s) / (ρ_g × v_th)
    ///
    /// # Stokes regime (s > 9λ/4, Re < 1)
    /// t_s = (2 ρ_m s²) / (9 μ)
    ///     = (2 ρ_m s²) / (9 × (1/3) ρ_g v_th λ)
    ///     = (2 ρ_m s²) / (3 ρ_g v_th λ)
    ///
    /// # References
    /// - Weidenschilling (1977), MNRAS 180, 57
    /// - Birnstiel et al. (2010), A&A 513, A79
    pub fn stopping_time<D: DiskModel>(&self, disk: &D, r: Length) -> Time {
        let rho_g = disk.midplane_density(r).to_grams_per_cm3();
        let v_th = disk.thermal_velocity(r).to_cm_per_sec();

        let t_s = match self.drag_regime(disk, r) {
            DragRegime::Epstein => {
                // t_s = ρ_m s / (ρ_g v_th)
                (self.material_density * self.size) / (rho_g * v_th)
            }
            DragRegime::Stokes => {
                // t_s = 2 ρ_m s² / (9 μ) where μ = (1/3) ρ_g v_th λ
                // Simplifies to: t_s = 2 ρ_m s² / (3 ρ_g v_th λ)
                let mfp = disk.mean_free_path(r).to_cm();
                (2.0 * self.material_density * self.size.powi(2)) / (3.0 * rho_g * v_th * mfp)
            }
            DragRegime::Turbulent => {
                // Not implemented - fall back to Stokes
                let mfp = disk.mean_free_path(r).to_cm();
                (2.0 * self.material_density * self.size.powi(2)) / (3.0 * rho_g * v_th * mfp)
            }
        };

        Time::from_seconds(t_s)
    }

    /// Stokes number: dimensionless stopping time.
    ///
    /// τ_s = t_s × Ω_K
    ///
    /// The Stokes number determines aerodynamic behavior:
    /// - τ_s << 1: particle tightly coupled to gas, moves with gas
    /// - τ_s ~ 1: maximum radial drift velocity
    /// - τ_s >> 1: particle decoupled from gas, Keplerian orbit
    ///
    /// For Epstein drag in the midplane:
    /// τ_s ≈ (π/2) × (ρ_m s) / Σ_g
    pub fn stokes_number<D: DiskModel>(&self, disk: &D, r: Length) -> f64 {
        let t_s = self.stopping_time(disk, r).to_seconds();
        let omega = disk.orbital_frequency(r).to_rad_per_sec();
        t_s * omega
    }

    /// Radial drift velocity due to gas drag.
    ///
    /// Particles feel a headwind from sub-Keplerian gas and drift inward.
    /// The drift velocity depends on the Stokes number:
    ///
    /// v_r = -2 η v_K τ_s / (1 + τ_s²)
    ///
    /// Peak drift occurs at τ_s = 1, where v_r = -η v_K.
    /// For typical disks, η ~ 0.003, so peak drift ~ 50 m/s.
    ///
    /// # References
    /// - Weidenschilling (1977), MNRAS 180, 57
    /// - Nakagawa et al. (1986), Icarus 67, 375
    pub fn radial_drift_velocity<D: DiskModel>(&self, disk: &D, r: Length) -> Velocity {
        let tau = self.stokes_number(disk, r);
        let eta = disk.pressure_gradient_parameter(r);
        let v_k = disk.keplerian_velocity(r).to_cm_per_sec();

        // v_r = -2 η v_K τ_s / (1 + τ_s²)
        let v_r = -2.0 * eta * v_k * tau / (1.0 + tau.powi(2));

        Velocity::from_cm_per_sec(v_r)
    }

    /// Azimuthal velocity deviation from Keplerian.
    ///
    /// Due to drag, particles orbit slightly slower than Keplerian:
    ///
    /// Δv_φ = v_K - v_φ = η v_K / (1 + τ_s²)
    ///
    /// Small particles (τ_s << 1) move with the sub-Keplerian gas.
    /// Large particles (τ_s >> 1) move at nearly Keplerian velocity.
    pub fn azimuthal_drift_velocity<D: DiskModel>(&self, disk: &D, r: Length) -> Velocity {
        let tau = self.stokes_number(disk, r);
        let eta = disk.pressure_gradient_parameter(r);
        let v_k = disk.keplerian_velocity(r).to_cm_per_sec();

        // Δv_φ = η v_K / (1 + τ_s²)
        let delta_v = eta * v_k / (1.0 + tau.powi(2));

        Velocity::from_cm_per_sec(delta_v)
    }

    /// Vertical settling velocity at height z above midplane.
    ///
    /// Particles settle toward the midplane on a timescale controlled
    /// by the stopping time:
    ///
    /// v_z = -Ω_K² z × t_s = -τ_s Ω_K z
    ///
    /// Valid for τ_s << 1 (well-coupled particles).
    pub fn settling_velocity<D: DiskModel>(&self, disk: &D, r: Length, z: Length) -> Velocity {
        let tau = self.stokes_number(disk, r);
        let omega = disk.orbital_frequency(r).to_rad_per_sec();
        let z_cm = z.to_cm();

        // v_z = -τ_s Ω z
        let v_z = -tau * omega * z_cm;

        Velocity::from_cm_per_sec(v_z)
    }

    /// Equilibrium scale height of a particle layer.
    ///
    /// Settling is balanced by turbulent diffusion, giving:
    ///
    /// h_p / h_g = √(α / (α + τ_s))
    ///
    /// where α is the turbulent viscosity parameter.
    ///
    /// - Small particles (τ_s << α): h_p ≈ h_g (well-mixed with gas)
    /// - Large particles (τ_s >> α): h_p ≈ h_g √(α/τ_s) (thin layer)
    pub fn equilibrium_scale_height<D: DiskModel>(&self, disk: &D, r: Length) -> Length {
        let tau = self.stokes_number(disk, r);
        let alpha = disk.alpha();
        let h_g = disk.scale_height(r).to_cm();

        // h_p / h_g = √(α / (α + τ_s))
        let h_p = h_g * (alpha / (alpha + tau)).sqrt();

        Length::from_cm(h_p)
    }
}
