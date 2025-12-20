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
    // Size-resolved properties
    // =========================================================================

    /// Stokes number for each size bin.
    ///
    /// Returns (size, Stokes_number) pairs for each bin in the distribution.
    pub fn stokes_numbers<D: DiskModel>(&self, disk: &D) -> Vec<(Length, f64)> {
        let r = self.radial_center();
        let rho_m = self.material_density();

        self.size_distribution
            .bins()
            .into_iter()
            .map(|(size, _mass)| {
                let particle = Particle::new(size, rho_m);
                (size, particle.stokes_number(disk, r))
            })
            .collect()
    }

    /// Radial drift velocity for each size bin.
    ///
    /// Returns (size, velocity) pairs. Negative velocity indicates inward drift.
    pub fn drift_velocities<D: DiskModel>(&self, disk: &D) -> Vec<(Length, Velocity)> {
        let r = self.radial_center();
        let rho_m = self.material_density();

        self.size_distribution
            .bins()
            .into_iter()
            .map(|(size, _mass)| {
                let particle = Particle::new(size, rho_m);
                (size, particle.radial_drift_velocity(disk, r))
            })
            .collect()
    }

    /// Relative velocity between particles of two sizes.
    ///
    /// This is the key quantity for collision rates. Includes contributions from:
    /// - Differential radial drift
    /// - Differential azimuthal drift
    /// - Turbulent relative velocities (Ormel & Cuzzi 2007)
    ///
    /// # Arguments
    /// * `disk` - Gas disk model
    /// * `s1` - First particle size
    /// * `s2` - Second particle size
    pub fn relative_velocity<D: DiskModel>(&self, disk: &D, s1: Length, s2: Length) -> Velocity {
        let r = self.radial_center();
        let rho_m = self.material_density();

        let p1 = Particle::new(s1, rho_m);
        let p2 = Particle::new(s2, rho_m);

        // Differential radial drift
        let v_r1 = p1.radial_drift_velocity(disk, r).to_cm_per_sec();
        let v_r2 = p2.radial_drift_velocity(disk, r).to_cm_per_sec();
        let dv_r = (v_r1 - v_r2).abs();

        // Differential azimuthal drift
        let v_phi1 = p1.azimuthal_drift_velocity(disk, r).to_cm_per_sec();
        let v_phi2 = p2.azimuthal_drift_velocity(disk, r).to_cm_per_sec();
        let dv_phi = (v_phi1 - v_phi2).abs();

        // Turbulent relative velocity
        let dv_turb = self.turbulent_relative_velocity(disk, s1, s2);

        // Quadrature sum (independent sources add in quadrature)
        let dv_total = (dv_r.powi(2) + dv_phi.powi(2) + dv_turb.powi(2)).sqrt();

        Velocity::from_cm_per_sec(dv_total)
    }

    /// Turbulent relative velocity between two particle sizes.
    ///
    /// Follows the prescription of Ormel & Cuzzi (2007) for closed-form
    /// turbulent collision velocities.
    ///
    /// For particles with Stokes numbers τ1 and τ2 in a turbulent flow:
    /// - If both τ << 1: particles couple to same eddies, low relative velocity
    /// - If τ1 ~ 1, τ2 << 1: large particle decoupled, small coupled → high Δv
    /// - If both τ >> 1: both decoupled from gas, low relative velocity
    fn turbulent_relative_velocity<D: DiskModel>(&self, disk: &D, s1: Length, s2: Length) -> f64 {
        let r = self.radial_center();
        let rho_m = self.material_density();
        let alpha = disk.alpha();

        let p1 = Particle::new(s1, rho_m);
        let p2 = Particle::new(s2, rho_m);

        let tau1 = p1.stokes_number(disk, r);
        let tau2 = p2.stokes_number(disk, r);

        // Turbulent velocity of the gas: v_turb ~ sqrt(α) × c_s
        let c_s = disk.sound_speed(r).to_cm_per_sec();
        let v_turb_gas = alpha.sqrt() * c_s;

        // Ormel & Cuzzi (2007) approximation for relative velocity
        // For simplicity, use the intermediate regime formula:
        // Δv_turb ≈ v_turb × sqrt(|τ1 - τ2| / (τ1 + τ2 + Re_t^(-1/2)))
        //
        // where Re_t is the turbulent Reynolds number.
        // For α ~ 10^-3 disks, Re_t ~ 10^8, so Re_t^(-1/2) ~ 10^-4 is negligible.
        //
        // Simplified formula:
        // Δv_turb ≈ v_turb × sqrt(3 × τ_s) for τ_s << 1
        // Δv_turb ≈ v_turb × sqrt(3 / τ_s) for τ_s >> 1
        //
        // General interpolation (Ormel & Cuzzi eq. 27):
        let tau_max = tau1.max(tau2);
        let tau_min = tau1.min(tau2);

        if tau_max < 1e-10 {
            // Both extremely small - Brownian motion dominates (not included here)
            return 0.0;
        }

        // Regime-dependent formula
        if tau_max < 1.0 {
            // Both tightly coupled: Δv ~ v_turb × sqrt(τ_max - τ_min)
            v_turb_gas * (1.5 * (tau_max - tau_min).abs()).sqrt()
        } else if tau_min > 1.0 {
            // Both loosely coupled: Δv ~ v_turb / sqrt(τ_min)
            v_turb_gas * (1.0 / tau_min).sqrt()
        } else {
            // Mixed regime: large particle decoupled, small coupled
            // This gives maximum relative velocity
            v_turb_gas * (1.5 * tau_max).sqrt().min(1.5 * v_turb_gas)
        }
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

    // =========================================================================
    // Gravitational Instability
    // =========================================================================

    /// Toomre Q parameter for the particle layer.
    ///
    /// The Toomre Q parameter measures gravitational stability:
    ///
    /// ```text
    /// Q = σ × Ω / (π × G × Σ)
    /// ```
    ///
    /// where:
    /// - σ is the velocity dispersion
    /// - Ω is the orbital frequency
    /// - G is the gravitational constant
    /// - Σ is the surface density
    ///
    /// **Stability criterion:**
    /// - Q > Q_crit (≈ 1-2): Gravitationally stable
    /// - Q < Q_crit: Gravitationally unstable, layer can fragment
    ///
    /// For a particle layer to collapse into planetesimals via the
    /// Goldreich-Ward mechanism or streaming instability, Q must be
    /// below the critical value.
    ///
    /// # References
    /// - Toomre (1964) - Original stability criterion for gaseous disks
    /// - Goldreich & Ward (1973) - Application to particle layers
    /// - Youdin & Goodman (2005) - Streaming instability
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::GasDisk;
    /// use stellar_forge::particles::ParticleBin;
    /// use units::Length;
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    /// let bin = ParticleBin::from_disk(&disk, r, width);
    ///
    /// let q = bin.toomre_q(&disk);
    /// // Initially well-mixed particles have high Q (stable)
    /// assert!(q > 10.0);
    /// ```
    pub fn toomre_q<D: DiskModel>(&self, disk: &D) -> f64 {
        use crate::disk::constants::G;

        let r = self.radial_center();
        let sigma = self.velocity_dispersion;
        let omega = disk.orbital_frequency(r).to_rad_per_sec();
        let surface_density = self.surface_density;

        // Q = σ × Ω / (π × G × Σ)
        sigma * omega / (PI * G * surface_density)
    }

    /// Richardson number for the particle layer.
    ///
    /// The Richardson number measures vertical shear stability:
    ///
    /// ```text
    /// Ri = (Ω_K² + N²) / (dv/dz)²
    /// ```
    ///
    /// where N² is the squared Brunt-Väisälä frequency (buoyancy oscillation).
    ///
    /// For a particle layer in a gas disk, the Richardson number must be
    /// sufficiently high (Ri > Ri_crit ≈ 0.25) to avoid Kelvin-Helmholtz
    /// instability from vertical shear.
    ///
    /// **Key physics:**
    /// - If particles settle too much, they create a sharp density gradient
    /// - Sharp gradients can trigger shear instabilities
    /// - These instabilities stir the layer, preventing further settling
    /// - Result: self-regulation of particle layer thickness
    ///
    /// This prevents unrealistic collapse in simulations and captures the
    /// physics of turbulence generation by particle layers.
    ///
    /// # Simplified Calculation
    ///
    /// For a thin particle layer in a gaseous disk, we use a simplified
    /// estimate based on the particle-to-gas density ratio and vertical
    /// structure:
    ///
    /// ```text
    /// Ri ≈ (h_p / h_g)² × (1 + ρ_p / ρ_g)
    /// ```
    ///
    /// This captures the essential physics:
    /// - Thicker particle layers (h_p → h_g) → high Ri → stable
    /// - Thin particle layers (h_p << h_g) → low Ri → potentially unstable
    /// - High particle density → stabilizing (more inertia)
    ///
    /// # References
    /// - Sekiya (1998) - Gravitational instability in particle disks
    /// - Youdin & Shu (2002) - Particle layer instabilities
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::GasDisk;
    /// use stellar_forge::particles::ParticleBin;
    /// use units::Length;
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    /// let bin = ParticleBin::from_disk(&disk, r, width);
    ///
    /// let ri = bin.richardson_number(&disk);
    /// // Well-mixed particles have high Ri (stable against shear)
    /// assert!(ri > 1.0);
    /// ```
    pub fn richardson_number<D: DiskModel>(&self, disk: &D) -> f64 {
        let r = self.radial_center();

        // Particle layer scale height
        let h_p = self.scale_height;

        // Gas disk scale height
        let h_g = disk.scale_height(r).to_cm();

        // Midplane densities
        let rho_p = self.midplane_density().to_grams_per_cm3();
        let rho_g = disk.midplane_density(r).to_grams_per_cm3();

        // Simplified Richardson number estimate
        // Ri ≈ (h_p / h_g)² × (1 + ρ_p / ρ_g)
        let height_ratio = h_p / h_g.max(1e-10);
        let density_enhancement = 1.0 + rho_p / rho_g.max(1e-10);

        height_ratio.powi(2) * density_enhancement
    }

    /// Check if the particle layer is gravitationally unstable.
    ///
    /// A particle layer can collapse into planetesimals if two conditions
    /// are met simultaneously:
    ///
    /// 1. **Gravitational instability**: Toomre Q < Q_crit
    ///    - Self-gravity overcomes random velocities
    ///    - Layer wants to clump
    ///
    /// 2. **Shear stability**: Richardson Ri > Ri_crit
    ///    - No self-excited turbulence from shear
    ///    - Layer remains coherent
    ///
    /// If both conditions hold, the layer is in the "Goldreich-Ward unstable"
    /// regime and can fragment into planetesimals.
    ///
    /// # Critical Values
    ///
    /// - **Q_crit ≈ 1.0-2.0**: Standard threshold from linear stability analysis
    /// - **Ri_crit ≈ 0.25**: Classic fluid dynamics threshold
    ///
    /// # Physical Interpretation
    ///
    /// - **Q < Q_crit, Ri > Ri_crit**: Goldreich-Ward unstable → planetesimals form
    /// - **Q < Q_crit, Ri < Ri_crit**: Self-excited turbulence → layer thickens, Q increases
    /// - **Q > Q_crit**: Gravitationally stable → no fragmentation
    ///
    /// This creates a self-regulating system where particle layers can only
    /// fragment in specific locations (e.g., pressure bumps, snow line) where
    /// conditions are just right.
    ///
    /// # References
    /// - Goldreich & Ward (1973) - Gravitational instability mechanism
    /// - Youdin & Goodman (2005) - Streaming instability
    /// - Johansen et al. (2007) - Numerical simulations
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::GasDisk;
    /// use stellar_forge::particles::ParticleBin;
    /// use units::Length;
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    /// let bin = ParticleBin::from_disk(&disk, r, width);
    ///
    /// // Initially stable (well-mixed, low density)
    /// assert!(!bin.is_gravitationally_unstable(&disk));
    /// ```
    pub fn is_gravitationally_unstable<D: DiskModel>(&self, disk: &D) -> bool {
        // Critical values from linear stability analysis
        const Q_CRIT: f64 = 1.5; // Conservative value between 1.0 and 2.0
        const RI_CRIT: f64 = 0.25; // Classic Richardson criterion

        let q = self.toomre_q(disk);
        let ri = self.richardson_number(disk);

        // Must satisfy both conditions:
        // 1. Gravitationally unstable (low Q)
        // 2. Shear stable (high Ri)
        q < Q_CRIT && ri > RI_CRIT
    }

    // =========================================================================
    // Planetesimal Formation
    // =========================================================================

    /// Jeans mass for the particle layer.
    ///
    /// When a particle layer becomes gravitationally unstable, it fragments
    /// into clumps with a characteristic mass scale set by the Jeans length.
    ///
    /// For a self-gravitating particle layer with surface density Σ and
    /// velocity dispersion σ:
    ///
    /// ```text
    /// λ_J = σ² / (G × Σ)
    /// M_J = π × λ_J² × Σ
    /// ```
    ///
    /// Simplifying:
    ///
    /// ```text
    /// M_J = π × σ⁴ / (G² × Σ)
    /// ```
    ///
    /// This is the characteristic mass of planetesimals that form from
    /// gravitational instability.
    ///
    /// # Physical Interpretation
    ///
    /// - **High velocity dispersion** → larger Jeans mass (harder to compress)
    /// - **High surface density** → smaller Jeans mass (stronger self-gravity)
    /// - Typical values: 10¹⁸ - 10²¹ g (1-100 km diameter)
    ///
    /// # References
    /// - Goldreich & Ward (1973) - Original calculation
    /// - Youdin & Shu (2002) - Particle layer physics
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::{DiskModel, GasDisk};
    /// use stellar_forge::particles::ParticleBin;
    /// use units::{Length, SurfaceDensity, Density};
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    ///
    /// // Create a dense, settled layer
    /// let mut bin = ParticleBin::from_disk(&disk, r, width);
    /// bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(100.0));
    /// bin.set_scale_height(disk.scale_height(r) * 0.1);
    /// bin.set_velocity_dispersion(disk.sound_speed(r) * 0.01);
    ///
    /// let m_j = bin.jeans_mass();
    ///
    /// // Typical planetesimal mass range
    /// assert!(m_j.to_grams() > 1e18);
    /// ```
    pub fn jeans_mass(&self) -> Mass {
        use crate::disk::constants::G;

        let sigma = self.velocity_dispersion;
        let surface_density = self.surface_density;

        // M_J = π × σ⁴ / (G² × Σ)
        let m_j = PI * sigma.powi(4) / (G.powi(2) * surface_density);

        Mass::from_grams(m_j)
    }

    /// Free-fall time for gravitational collapse.
    ///
    /// When a particle layer becomes unstable, it collapses on the free-fall
    /// timescale:
    ///
    /// ```text
    /// t_ff = 1 / sqrt(G × ρ)
    /// ```
    ///
    /// where ρ is the midplane volume density of the particle layer.
    ///
    /// This sets the timescale for planetesimal formation once instability
    /// criteria are met.
    ///
    /// # Physical Interpretation
    ///
    /// - Denser layers collapse faster
    /// - Typical values: 10-100 orbits
    /// - Much faster than coagulation timescales
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::{DiskModel, GasDisk};
    /// use stellar_forge::particles::ParticleBin;
    /// use units::{Length, SurfaceDensity};
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    ///
    /// let mut bin = ParticleBin::from_disk(&disk, r, width);
    /// bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(100.0));
    /// bin.set_scale_height(disk.scale_height(r) * 0.1);
    ///
    /// let t_ff = bin.freefall_time();
    /// let t_orb = disk.orbital_period(r);
    ///
    /// // Free-fall time should be comparable to orbital period
    /// let ratio = t_ff.to_seconds() / t_orb.to_seconds();
    /// assert!(ratio > 0.1 && ratio < 100.0);
    /// ```
    pub fn freefall_time(&self) -> Time {
        use crate::disk::constants::G;

        let rho = self.midplane_density().to_grams_per_cm3();

        // t_ff = sqrt(3π / (32 G ρ))
        // Simplified to 1 / sqrt(G × ρ) for order-of-magnitude
        let t_ff = 1.0 / (G * rho).sqrt();

        Time::from_seconds(t_ff)
    }

    /// Characteristic planetesimal size from Jeans mass.
    ///
    /// Converts the Jeans mass to a physical radius assuming spherical
    /// geometry and material density.
    ///
    /// ```text
    /// R = (3 M_J / (4π ρ_m))^(1/3)
    /// ```
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::{DiskModel, GasDisk};
    /// use stellar_forge::particles::ParticleBin;
    /// use units::{Length, SurfaceDensity};
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    ///
    /// let mut bin = ParticleBin::from_disk(&disk, r, width);
    /// bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(100.0));
    /// bin.set_scale_height(disk.scale_height(r) * 0.1);
    /// bin.set_velocity_dispersion(disk.sound_speed(r) * 0.01);
    ///
    /// let size = bin.planetesimal_size();
    ///
    /// // Should be km-scale
    /// assert!(size.to_km() > 1.0);
    /// ```
    pub fn planetesimal_size(&self) -> Length {
        let m_j = self.jeans_mass().to_grams();
        let rho_m = self.material_density().to_grams_per_cm3();

        // R = (3 M / (4π ρ))^(1/3)
        let radius = ((3.0 * m_j) / (4.0 * PI * rho_m)).powf(1.0 / 3.0);

        Length::from_cm(radius)
    }

    /// Attempt to form planetesimals from gravitational instability.
    ///
    /// Checks if the particle bin is gravitationally unstable and, if so,
    /// converts a fraction of its mass into planetesimals.
    ///
    /// # Process
    ///
    /// 1. Check stability criteria (Toomre Q and Richardson number)
    /// 2. If unstable, calculate Jeans mass and formation efficiency
    /// 3. Remove mass from the particle bin
    /// 4. Return a formation event describing the planetesimals created
    ///
    /// # Arguments
    /// * `disk` - Gas disk model for stability checks
    /// * `rng` - Random number generator for stochastic efficiency
    ///
    /// # Returns
    /// - `Some(event)` if planetesimals formed
    /// - `None` if bin is stable
    ///
    /// # Example
    /// ```
    /// use stellar_forge::disk::{DiskModel, GasDisk};
    /// use stellar_forge::particles::ParticleBin;
    /// use units::{Length, SurfaceDensity};
    /// use rand::SeedableRng;
    /// use rand_chacha::ChaChaRng;
    ///
    /// let disk = GasDisk::mmsn();
    /// let r = Length::from_au(5.0);
    /// let width = Length::from_au(0.5);
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    ///
    /// // Create an unstable bin
    /// let mut bin = ParticleBin::from_disk(&disk, r, width);
    /// bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(100.0));
    /// bin.set_scale_height(disk.scale_height(r) * 0.5);
    /// bin.set_velocity_dispersion(disk.sound_speed(r) * 0.001);
    ///
    /// // Attempt formation
    /// if let Some(event) = bin.attempt_planetesimal_formation(&disk, &mut rng) {
    ///     println!("Formed {} planetesimals at {} AU",
    ///              event.number_formed(),
    ///              event.location.to_au());
    /// }
    /// ```
    pub fn attempt_planetesimal_formation<D: DiskModel>(
        &mut self,
        disk: &D,
        rng: &mut rand_chacha::ChaChaRng,
    ) -> Option<crate::particles::PlanetesimalFormationEvent> {
        use crate::particles::PlanetesimalFormationEvent;

        // Check if gravitationally unstable
        if !self.is_gravitationally_unstable(disk) {
            return None;
        }

        // Calculate properties for planetesimal formation
        let jeans_mass = self.jeans_mass();
        let jeans_size = self.planetesimal_size();
        let available_mass = self.total_mass();
        let freefall_time = self.freefall_time();

        // Sample formation efficiency based on how unstable the layer is
        let q = self.toomre_q(disk);
        let efficiency = PlanetesimalFormationEvent::sample_efficiency(q, rng);

        // Create formation event
        let event = PlanetesimalFormationEvent::new(
            self.radial_center(),
            jeans_mass,
            jeans_size,
            available_mass,
            freefall_time,
            efficiency,
            self.material_density(),
        );

        // Remove mass from particle bin (it's now in planetesimals)
        let mass_converted = event.total_mass.to_grams();
        let area = self.area();
        let new_surface_density = self.surface_density - mass_converted / area;

        // Ensure non-negative surface density
        self.surface_density = new_surface_density.max(0.0);

        Some(event)
    }
}
