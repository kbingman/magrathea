//! Grid-based protoplanetary disk model.
//!
//! Stores surface density on a discrete radial grid, enabling:
//! - Viscous evolution via diffusion equation
//! - Local modifications (gap opening, accretion)
//! - Arbitrary density profiles
//!
//! The grid uses logarithmic spacing for uniform resolution in log(r).

use units::{Length, Mass, SurfaceDensity, Temperature};

use super::disk_model::{DiskMass, DiskModel};
use super::gas_disk::GasDisk;
use crate::disk::constants::PI;

/// A protoplanetary disk with surface density stored on a radial grid.
///
/// Uses logarithmic spacing for the radial grid and log-linear interpolation
/// for queries between grid points.
#[derive(Debug, Clone)]
pub struct GridDisk {
    /// Radial grid points (cm), logarithmically spaced
    radii: Vec<f64>,

    /// Surface density at each grid point (g/cm²)
    sigma: Vec<f64>,

    /// Temperature at reference radius (K)
    temp_0: f64,

    /// Temperature power-law exponent (T ∝ r^(-q))
    temp_exponent: f64,

    /// Reference radius for temperature profile (cm)
    r_0: f64,

    /// Central stellar mass (g)
    stellar_mass: f64,

    /// Shakura-Sunyaev viscosity parameter
    alpha: f64,
}

impl GridDisk {
    /// Create a new grid disk with the given parameters.
    ///
    /// # Arguments
    /// * `radii` - Radial grid points (should be logarithmically spaced)
    /// * `sigma` - Surface density at each grid point
    /// * `temp_0` - Temperature at reference radius
    /// * `temp_exponent` - Temperature power-law exponent
    /// * `r_0` - Reference radius for temperature
    /// * `stellar_mass` - Central stellar mass
    /// * `alpha` - Viscosity parameter
    pub fn new(
        radii: Vec<f64>,
        sigma: Vec<f64>,
        temp_0: f64,
        temp_exponent: f64,
        r_0: f64,
        stellar_mass: f64,
        alpha: f64,
    ) -> Self {
        assert_eq!(
            radii.len(),
            sigma.len(),
            "radii and sigma must have same length"
        );
        assert!(radii.len() >= 2, "need at least 2 grid points");

        Self {
            radii,
            sigma,
            temp_0,
            temp_exponent,
            r_0,
            stellar_mass,
            alpha,
        }
    }

    /// Create a grid disk from a power-law profile.
    ///
    /// # Arguments
    /// * `stellar_mass` - Central stellar mass
    /// * `sigma_0` - Surface density at r_0 (g/cm²)
    /// * `sigma_exponent` - Surface density power-law exponent
    /// * `temp_0` - Temperature at r_0 (K)
    /// * `temp_exponent` - Temperature power-law exponent
    /// * `r_0` - Reference radius
    /// * `inner_radius` - Inner edge of disk
    /// * `outer_radius` - Outer edge of disk
    /// * `n_radii` - Number of grid points
    /// * `alpha` - Viscosity parameter
    #[allow(clippy::too_many_arguments)]
    pub fn from_power_law(
        stellar_mass: Mass,
        sigma_0: f64,
        sigma_exponent: f64,
        temp_0: f64,
        temp_exponent: f64,
        r_0: Length,
        inner_radius: Length,
        outer_radius: Length,
        n_radii: usize,
        alpha: f64,
    ) -> Self {
        let radii = log_spaced_grid(inner_radius.to_cm(), outer_radius.to_cm(), n_radii);
        let r_0_cm = r_0.to_cm();

        let sigma: Vec<f64> = radii
            .iter()
            .map(|&r| sigma_0 * (r / r_0_cm).powf(-sigma_exponent))
            .collect();

        Self {
            radii,
            sigma,
            temp_0,
            temp_exponent,
            r_0: r_0_cm,
            stellar_mass: stellar_mass.to_grams(),
            alpha,
        }
    }

    /// Create a grid disk from an existing GasDisk.
    ///
    /// This samples the power-law disk onto the grid.
    pub fn from_gas_disk(disk: &GasDisk, n_radii: usize) -> Self {
        let radii = log_spaced_grid(
            disk.inner_radius.to_cm(),
            disk.outer_radius.to_cm(),
            n_radii,
        );

        let sigma: Vec<f64> = radii
            .iter()
            .map(|&r| disk.surface_density(Length::from_cm(r)).to_grams_per_cm2())
            .collect();

        Self {
            radii,
            sigma,
            temp_0: disk.temperature_0.to_kelvin(),
            temp_exponent: disk.temp_exponent,
            r_0: disk.r_0.to_cm(),
            stellar_mass: disk.stellar_mass.to_grams(),
            alpha: disk.alpha,
        }
    }

    /// Number of grid points.
    pub fn n_radii(&self) -> usize {
        self.radii.len()
    }

    /// Get the radial grid points.
    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    /// Get the surface density values.
    pub fn sigma(&self) -> &[f64] {
        &self.sigma
    }

    /// Get mutable access to surface density values.
    ///
    /// Use this for viscous evolution or other modifications.
    pub fn sigma_mut(&mut self) -> &mut [f64] {
        &mut self.sigma
    }

    /// Interpolate surface density at radius r using log-linear interpolation.
    ///
    /// Linear interpolation in log-log space matches power-law profiles exactly.
    fn interpolate_sigma(&self, r: f64) -> f64 {
        let log_r = r.ln();

        // Find bracketing indices
        let i = self
            .radii
            .iter()
            .position(|&x| x.ln() >= log_r)
            .unwrap_or(self.radii.len())
            .saturating_sub(1)
            .min(self.radii.len() - 2);

        let log_r_i = self.radii[i].ln();
        let log_r_ip1 = self.radii[i + 1].ln();

        // Linear interpolation parameter
        let t = (log_r - log_r_i) / (log_r_ip1 - log_r_i);

        // Interpolate in log space
        let log_sigma_i = self.sigma[i].ln();
        let log_sigma_ip1 = self.sigma[i + 1].ln();

        (log_sigma_i + t * (log_sigma_ip1 - log_sigma_i)).exp()
    }
}

// =============================================================================
// DiskModel trait implementation
// =============================================================================

impl DiskModel for GridDisk {
    fn surface_density(&self, r: Length) -> SurfaceDensity {
        SurfaceDensity::from_grams_per_cm2(self.interpolate_sigma(r.to_cm()))
    }

    fn temperature(&self, r: Length) -> Temperature {
        let ratio = r.to_cm() / self.r_0;
        Temperature::from_kelvin(self.temp_0 * ratio.powf(-self.temp_exponent))
    }

    fn stellar_mass(&self) -> Mass {
        Mass::from_grams(self.stellar_mass)
    }

    fn alpha(&self) -> f64 {
        self.alpha
    }

    fn inner_radius(&self) -> Length {
        Length::from_cm(self.radii[0])
    }

    fn outer_radius(&self) -> Length {
        Length::from_cm(*self.radii.last().unwrap())
    }

    // Uses default numerical pressure_gradient_log
}

impl DiskMass for GridDisk {
    /// Total disk mass using trapezoidal integration.
    fn total_mass(&self) -> Mass {
        let mut mass = 0.0;

        for i in 0..self.radii.len() - 1 {
            let r1 = self.radii[i];
            let r2 = self.radii[i + 1];
            let s1 = self.sigma[i];
            let s2 = self.sigma[i + 1];

            // Trapezoidal rule: ∫ 2πr Σ dr ≈ π(r2² - r1²) × (Σ1 + Σ2) / 2
            // More accurate for power-law: use geometric mean
            let r_mid = (r1 * r2).sqrt();
            let s_mid = (s1 * s2).sqrt();
            let dr = r2 - r1;

            mass += 2.0 * PI * r_mid * s_mid * dr;
        }

        Mass::from_grams(mass)
    }
}

// =============================================================================
// Helper functions
// =============================================================================

/// Generate a logarithmically spaced grid.
fn log_spaced_grid(r_min: f64, r_max: f64, n: usize) -> Vec<f64> {
    let log_min = r_min.ln();
    let log_max = r_max.ln();

    (0..n)
        .map(|i| {
            let frac = i as f64 / (n - 1) as f64;
            (log_min + frac * (log_max - log_min)).exp()
        })
        .collect()
}
