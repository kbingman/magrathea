//! Spatial bins for particle populations.
//!
//! A `ParticleBin` represents a population of solid particles in a radial
//! annulus of the disk. It combines:
//! - Spatial location (radial position and width)
//! - Size distribution (how mass is distributed across particle sizes)
//! - Bulk properties (surface density, scale height, velocity dispersion)
//!
//! This enables efficient statistical treatment of ~10^15 particles without
//! tracking individuals.

use units::{Density, Length, Mass, SurfaceDensity, Time, Velocity};

use crate::disk::DiskModel;
use crate::disk::constants::PI;
use crate::particles::{Particle, SizeDistribution};

/// Default dust-to-gas mass ratio (solar metallicity).
const DUST_TO_GAS_RATIO: f64 = 0.01;

/// Default material density for silicate grains (g/cm³).
const SILICATE_DENSITY: f64 = 3.0;

/// Number of size bins for discretized distributions.
const N_SIZE_BINS: usize = 30;

/// A population of solid particles in a radial annulus of the disk.
///
/// Tracks the spatial distribution, size distribution, and bulk properties
/// of particles at a given location in the disk.
#[derive(Debug, Clone)]
pub struct ParticleBin {
    /// Center of the radial annulus (cm)
    radial_center: f64,

    /// Width of the radial annulus (cm)
    radial_width: f64,

    /// Size distribution of particles
    size_distribution: SizeDistribution,

    /// Surface density of solids (g/cm²)
    surface_density: f64,

    /// Vertical scale height of particle layer (cm)
    scale_height: f64,

    /// Random velocity dispersion (cm/s)
    velocity_dispersion: f64,
}

impl ParticleBin {
    /// Create a particle bin with explicit parameters.
    ///
    /// # Arguments
    /// * `radial_center` - Center of the radial annulus
    /// * `radial_width` - Width of the annulus
    /// * `size_distribution` - Distribution of particle sizes
    /// * `surface_density` - Total surface density of solids
    /// * `scale_height` - Vertical scale height
    pub fn new(
        radial_center: Length,
        radial_width: Length,
        size_distribution: SizeDistribution,
        surface_density: SurfaceDensity,
        scale_height: Length,
    ) -> Self {
        Self {
            radial_center: radial_center.to_cm(),
            radial_width: radial_width.to_cm(),
            size_distribution,
            surface_density: surface_density.to_grams_per_cm2(),
            scale_height: scale_height.to_cm(),
            velocity_dispersion: 0.0,
        }
    }

    /// Create a particle bin from disk conditions at radius r.
    ///
    /// Uses standard assumptions:
    /// - Dust-to-gas ratio: 0.01 (solar metallicity)
    /// - Initial size distribution: MRN (0.1 μm to 1 μm, q = 3.5)
    /// - Material density: 3.0 g/cm³ (silicate)
    /// - Initial scale height: same as gas (well-mixed)
    /// - Velocity dispersion: 1% of sound speed
    ///
    /// # Arguments
    /// * `disk` - The gas disk model
    /// * `r` - Radial position
    /// * `width` - Width of the annulus
    pub fn from_disk<D: DiskModel>(disk: &D, r: Length, width: Length) -> Self {
        let r_cm = r.to_cm();
        let width_cm = width.to_cm();

        // Surface density from dust-to-gas ratio
        let sigma_gas = disk.surface_density(r).to_grams_per_cm2();
        let sigma_dust = sigma_gas * DUST_TO_GAS_RATIO;

        // Calculate annulus area for total mass
        let r_in = r_cm - width_cm / 2.0;
        let r_out = r_cm + width_cm / 2.0;
        let area = PI * (r_out.powi(2) - r_in.powi(2));
        let total_mass = sigma_dust * area;

        // MRN size distribution, converted to binned for uniform iteration
        let size_dist = SizeDistribution::mrn(
            Mass::from_grams(total_mass),
            Density::from_grams_per_cm3(SILICATE_DENSITY),
        )
        .to_binned(N_SIZE_BINS);

        // Initially well-mixed with gas
        let h_gas = disk.scale_height(r).to_cm();

        // Small random velocities (fraction of sound speed)
        let c_s = disk.sound_speed(r).to_cm_per_sec();
        let v_disp = 0.01 * c_s;

        Self {
            radial_center: r_cm,
            radial_width: width_cm,
            size_distribution: size_dist,
            surface_density: sigma_dust,
            scale_height: h_gas,
            velocity_dispersion: v_disp,
        }
    }

    // =========================================================================
    // Geometric queries
    // =========================================================================

    /// Center of the radial annulus.
    pub fn radial_center(&self) -> Length {
        Length::from_cm(self.radial_center)
    }

    /// Width of the radial annulus.
    pub fn radial_width(&self) -> Length {
        Length::from_cm(self.radial_width)
    }

    /// Inner edge of the annulus.
    pub fn inner_radius(&self) -> Length {
        Length::from_cm(self.radial_center - self.radial_width / 2.0)
    }

    /// Outer edge of the annulus.
    pub fn outer_radius(&self) -> Length {
        Length::from_cm(self.radial_center + self.radial_width / 2.0)
    }

    /// Area of the annulus: π(r_out² - r_in²)
    pub fn area(&self) -> f64 {
        let r_in = self.radial_center - self.radial_width / 2.0;
        let r_out = self.radial_center + self.radial_width / 2.0;
        PI * (r_out.powi(2) - r_in.powi(2))
    }

    // =========================================================================
    // Mass and density queries
    // =========================================================================

    /// Surface density of solids.
    pub fn surface_density(&self) -> SurfaceDensity {
        SurfaceDensity::from_grams_per_cm2(self.surface_density)
    }

    /// Total mass in this bin: Σ × Area
    pub fn total_mass(&self) -> Mass {
        Mass::from_grams(self.surface_density * self.area())
    }

    /// Midplane volume density: Σ / (√(2π) × h)
    ///
    /// Assumes a Gaussian vertical density profile.
    pub fn midplane_density(&self) -> Density {
        let rho = self.surface_density / ((2.0 * PI).sqrt() * self.scale_height);
        Density::from_grams_per_cm3(rho)
    }

    /// Vertical scale height of the particle layer.
    pub fn scale_height(&self) -> Length {
        Length::from_cm(self.scale_height)
    }

    /// Random velocity dispersion.
    pub fn velocity_dispersion(&self) -> Velocity {
        Velocity::from_cm_per_sec(self.velocity_dispersion)
    }

    // =========================================================================
    // Size distribution access
    // =========================================================================

    /// Reference to the size distribution.
    pub fn size_distribution(&self) -> &SizeDistribution {
        &self.size_distribution
    }

    /// Mutable reference to the size distribution.
    pub fn size_distribution_mut(&mut self) -> &mut SizeDistribution {
        &mut self.size_distribution
    }

    /// Material density of particles.
    pub fn material_density(&self) -> Density {
        self.size_distribution.material_density()
    }

    /// Mean particle size (mass-weighted).
    pub fn mean_size(&self) -> Length {
        self.size_distribution.mean_size()
    }

    // =========================================================================
    // Aerodynamic properties
    // =========================================================================

    /// Mean stopping time (mass-weighted over size distribution).
    ///
    /// Computed by sampling the size distribution and averaging.
    pub fn mean_stopping_time<D: DiskModel>(&self, disk: &D) -> Time {
        let r = self.radial_center();
        let rho_m = self.material_density();

        // For now, use the mean size as representative
        // A more accurate implementation would integrate over the distribution
        let mean_s = self.mean_size();
        let particle = Particle::new(mean_s, rho_m);

        particle.stopping_time(disk, r)
    }

    /// Mean Stokes number (mass-weighted).
    ///
    /// τ_s = t_s × Ω_K
    pub fn mean_stokes_number<D: DiskModel>(&self, disk: &D) -> f64 {
        let r = self.radial_center();
        let rho_m = self.material_density();

        let mean_s = self.mean_size();
        let particle = Particle::new(mean_s, rho_m);

        particle.stokes_number(disk, r)
    }

    /// Mean radial drift velocity (mass-weighted).
    ///
    /// Negative values indicate inward drift.
    pub fn mean_drift_velocity<D: DiskModel>(&self, disk: &D) -> Velocity {
        let r = self.radial_center();
        let rho_m = self.material_density();

        let mean_s = self.mean_size();
        let particle = Particle::new(mean_s, rho_m);

        particle.radial_drift_velocity(disk, r)
    }

    /// Equilibrium scale height for the mean particle size.
    ///
    /// Computed from the balance of settling and turbulent diffusion.
    pub fn equilibrium_scale_height<D: DiskModel>(&self, disk: &D) -> Length {
        let r = self.radial_center();
        let rho_m = self.material_density();

        let mean_s = self.mean_size();
        let particle = Particle::new(mean_s, rho_m);

        particle.equilibrium_scale_height(disk, r)
    }

    // =========================================================================
    // Mutation
    // =========================================================================

    /// Set the surface density.
    pub fn set_surface_density(&mut self, sigma: SurfaceDensity) {
        self.surface_density = sigma.to_grams_per_cm2();
    }

    /// Set the scale height.
    pub fn set_scale_height(&mut self, h: Length) {
        self.scale_height = h.to_cm();
    }

    /// Set the velocity dispersion.
    pub fn set_velocity_dispersion(&mut self, v: Velocity) {
        self.velocity_dispersion = v.to_cm_per_sec();
    }

    /// Relax scale height toward equilibrium over timestep dt.
    ///
    /// Uses exponential relaxation: h → h_eq + (h - h_eq) × exp(-dt/t_settle)
    /// where t_settle ≈ 1 / (τ_s × Ω_K).
    pub fn relax_scale_height<D: DiskModel>(&mut self, disk: &D, dt: Time) {
        let r = self.radial_center();
        let h_eq = self.equilibrium_scale_height(disk).to_cm();
        let tau = self.mean_stokes_number(disk);
        let omega = disk.orbital_frequency(r).to_rad_per_sec();

        // Settling timescale: t_settle ≈ 1 / (τ_s × Ω_K)
        // For small τ_s, settling is slow; for large τ_s, settling is fast
        let t_settle = 1.0 / (tau.max(1e-10) * omega);

        // Exponential relaxation
        let dt_s = dt.to_seconds();
        let factor = (-dt_s / t_settle).exp();

        self.scale_height = h_eq + (self.scale_height - h_eq) * factor;
    }
}
