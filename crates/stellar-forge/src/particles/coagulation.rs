//! Particle coagulation (sticking) for dust growth.
//!
//! Models the growth of dust particles through collisions. When particles
//! collide at low velocities, they stick together forming larger aggregates.
//!
//! # Physics
//!
//! The evolution of the size distribution follows the Smoluchowski equation:
//!
//! ```text
//! dm_k/dt = (1/2) Σ_{i+j→k} K_ij × n_i × n_j × m_k
//!         - m_k × Σ_j K_kj × n_j
//! ```
//!
//! where K_ij is the collision kernel between size bins i and j.
//!
//! The collision kernel has the form:
//!
//! ```text
//! K_ij = π(s_i + s_j)² × Δv_ij
//! ```
//!
//! where Δv_ij is the relative velocity between particles of sizes s_i and s_j.
//!
//! # References
//!
//! - Smoluchowski (1916) - Original coagulation equation
//! - Brauer et al. (2008) - Dust coagulation in protoplanetary disks
//! - Birnstiel et al. (2010) - Gas and dust evolution in disks

use std::f64::consts::PI;

use units::{Length, Time};

use crate::disk::DiskModel;
use crate::particles::{ParticleBin, SizeDistribution};

/// Collision kernel calculator for particle coagulation.
///
/// Computes collision rates between size bins based on geometric cross-sections
/// and relative velocities.
pub struct Coagulation {
    /// Precomputed collision kernels K_ij (cm³/s per particle pair)
    kernel: Vec<Vec<f64>>,

    /// Number of bins
    n_bins: usize,

    /// Bin center sizes (cm)
    bin_sizes: Vec<f64>,

    /// Bin edges (cm), length n_bins + 1
    bin_edges: Vec<f64>,
}

impl Coagulation {
    /// Create a coagulation calculator for a particle bin.
    ///
    /// Precomputes the collision kernel matrix K_ij for all size bin pairs.
    ///
    /// # Arguments
    /// * `particle_bin` - The particle population (must have binned size distribution)
    /// * `disk` - Gas disk model for computing relative velocities
    pub fn new<D: DiskModel>(particle_bin: &ParticleBin, disk: &D) -> Self {
        let size_dist = particle_bin.size_distribution();
        let n_bins = size_dist.n_bins();

        assert!(n_bins > 0, "Size distribution must be binned");

        // Extract bin information
        let bins = size_dist.bins();
        let bin_sizes: Vec<f64> = bins.iter().map(|(s, _)| s.to_cm()).collect();

        // Get bin edges from distribution
        let bin_edges = Self::extract_bin_edges(size_dist);

        // Compute collision kernel matrix
        let mut kernel = vec![vec![0.0; n_bins]; n_bins];

        for i in 0..n_bins {
            for j in i..n_bins {
                let s_i = Length::from_cm(bin_sizes[i]);
                let s_j = Length::from_cm(bin_sizes[j]);

                // Geometric cross-section: π(s_i + s_j)²
                let cross_section = PI * (bin_sizes[i] + bin_sizes[j]).powi(2);

                // Relative velocity from ParticleBin
                let delta_v = particle_bin
                    .relative_velocity(disk, s_i, s_j)
                    .to_cm_per_sec();

                // Kernel: K_ij = σ × Δv
                let k_ij = cross_section * delta_v;

                kernel[i][j] = k_ij;
                kernel[j][i] = k_ij; // Symmetric
            }
        }

        Self {
            kernel,
            n_bins,
            bin_sizes,
            bin_edges,
        }
    }

    /// Extract bin edges from a size distribution.
    fn extract_bin_edges(dist: &SizeDistribution) -> Vec<f64> {
        match dist {
            SizeDistribution::Binned { bin_edges, .. } => bin_edges.clone(),
            SizeDistribution::Monodisperse { size, .. } => {
                vec![*size * 0.5, *size * 2.0]
            }
            SizeDistribution::PowerLaw { s_min, s_max, .. } => {
                vec![*s_min, *s_max]
            }
        }
    }

    /// Get the collision kernel between bins i and j.
    ///
    /// Returns K_ij in units of cm³/s (per particle pair).
    pub fn kernel(&self, i: usize, j: usize) -> f64 {
        self.kernel[i][j]
    }

    /// Find which bin a particle of given size falls into.
    ///
    /// Returns the bin index, or the last bin if size exceeds maximum.
    fn find_bin(&self, size_cm: f64) -> usize {
        for i in 0..self.n_bins {
            if size_cm < self.bin_edges[i + 1] {
                return i;
            }
        }
        self.n_bins - 1
    }

    /// Compute the size of a particle formed by two particles sticking.
    ///
    /// Mass is conserved: m_new = m_i + m_j
    /// For spheres: s_new = (s_i³ + s_j³)^(1/3)
    fn combined_size(&self, s_i: f64, s_j: f64) -> f64 {
        (s_i.powi(3) + s_j.powi(3)).powf(1.0 / 3.0)
    }

    /// Evolve the size distribution forward in time via sticking.
    ///
    /// Implements the Smoluchowski coagulation equation:
    /// - Pairs of particles collide and stick
    /// - Resulting particle goes into appropriate size bin
    /// - Total mass is conserved
    ///
    /// # Arguments
    /// * `dist` - Size distribution to evolve (must be Binned)
    /// * `surface_density` - Surface density of particles (g/cm²)
    /// * `scale_height` - Vertical scale height of particle layer (cm)
    /// * `dt` - Time step
    ///
    /// # Returns
    /// New size distribution after sticking, with conserved total mass.
    pub fn evolve_sticking(
        &self,
        dist: &SizeDistribution,
        surface_density: f64,
        scale_height: f64,
        dt: Time,
    ) -> SizeDistribution {
        let dt_s = dt.to_seconds();

        // Extract current state
        let (bin_edges, mass_per_bin, material_density) = match dist {
            SizeDistribution::Binned {
                bin_edges,
                mass_per_bin,
                material_density,
            } => (bin_edges.clone(), mass_per_bin.clone(), *material_density),
            _ => panic!("evolve_sticking requires a Binned distribution"),
        };

        let n = self.n_bins;

        // Total mass and mass fractions
        let total_mass: f64 = mass_per_bin.iter().sum();
        if total_mass == 0.0 {
            return dist.clone();
        }

        let mass_fraction: Vec<f64> = mass_per_bin.iter().map(|m| m / total_mass).collect();

        // Number density per unit volume at midplane: n_i = Σ_i / (sqrt(2π) H m_i)
        // where Σ_i = Σ_total × f_i
        let sqrt_2pi = (2.0 * PI).sqrt();
        let number_density: Vec<f64> = (0..n)
            .map(|i| {
                let s = self.bin_sizes[i];
                let m_particle = (4.0 / 3.0) * PI * material_density * s.powi(3);
                // Surface density in this bin
                let sigma_i = surface_density * mass_fraction[i];
                // Volume density
                sigma_i / (sqrt_2pi * scale_height * m_particle)
            })
            .collect();

        // Work in mass fraction space to avoid precision issues
        // df_k/dt = gains - losses, all normalized by total surface density
        let mut fraction_change = vec![0.0; n];

        for i in 0..n {
            for j in i..n {
                let n_i = number_density[i];
                let n_j = number_density[j];

                if n_i == 0.0 || n_j == 0.0 {
                    continue;
                }

                // Collision rate per unit volume: R_ij = K_ij × n_i × n_j
                // For i == j, we need factor of 1/2 (avoid double counting)
                let symmetry_factor = if i == j { 0.5 } else { 1.0 };
                let collision_rate = symmetry_factor * self.kernel[i][j] * n_i * n_j;

                // Mass involved per collision
                let s_i = self.bin_sizes[i];
                let s_j = self.bin_sizes[j];
                let m_i = (4.0 / 3.0) * PI * material_density * s_i.powi(3);
                let m_j = (4.0 / 3.0) * PI * material_density * s_j.powi(3);
                let m_combined = m_i + m_j;

                // Size of combined particle
                let s_combined = self.combined_size(s_i, s_j);
                let k = self.find_bin(s_combined);

                // Collision rate times dt gives number of collisions per unit volume
                let collisions_per_vol = collision_rate * dt_s;

                // Mass transferred per unit volume (g/cm³)
                // Each collision transfers m_i from bin i, m_j from bin j, m_combined to bin k
                let mass_loss_i_per_vol = collisions_per_vol * m_i;
                let mass_loss_j_per_vol = collisions_per_vol * m_j;
                let mass_gain_k_per_vol = collisions_per_vol * m_combined;

                // Convert to surface density change: dΣ = dρ × H × sqrt(2π)
                let loss_sigma_i = mass_loss_i_per_vol * scale_height * sqrt_2pi;
                let loss_sigma_j = mass_loss_j_per_vol * scale_height * sqrt_2pi;
                let gain_sigma_k = mass_gain_k_per_vol * scale_height * sqrt_2pi;

                // Convert to mass fraction change: df = dΣ / Σ_total
                let df_loss_i = loss_sigma_i / surface_density;
                let df_loss_j = loss_sigma_j / surface_density;
                let df_gain_k = gain_sigma_k / surface_density;

                // Limit to available mass fraction (stability)
                let max_df_i = 0.5 * mass_fraction[i];
                let max_df_j = 0.5 * mass_fraction[j];

                let scale_factor = if df_loss_i > 0.0 || df_loss_j > 0.0 {
                    let scale_i = if df_loss_i > max_df_i {
                        max_df_i / df_loss_i
                    } else {
                        1.0
                    };
                    let scale_j = if df_loss_j > max_df_j {
                        max_df_j / df_loss_j
                    } else {
                        1.0
                    };
                    scale_i.min(scale_j)
                } else {
                    1.0
                };

                // Apply scaled changes
                fraction_change[i] -= df_loss_i * scale_factor;
                if j != i {
                    fraction_change[j] -= df_loss_j * scale_factor;
                }
                fraction_change[k] += df_gain_k * scale_factor;
            }
        }

        // Apply fraction changes
        let mut new_fraction: Vec<f64> = mass_fraction
            .iter()
            .zip(fraction_change.iter())
            .map(|(&f, &df)| (f + df).max(0.0))
            .collect();

        // Renormalize to sum to 1
        let frac_sum: f64 = new_fraction.iter().sum();
        if frac_sum > 0.0 {
            for f in &mut new_fraction {
                *f /= frac_sum;
            }
        }

        // Convert back to absolute mass
        let new_mass: Vec<f64> = new_fraction.iter().map(|&f| f * total_mass).collect();

        SizeDistribution::Binned {
            bin_edges,
            mass_per_bin: new_mass,
            material_density,
        }
    }

    /// Number of size bins.
    pub fn n_bins(&self) -> usize {
        self.n_bins
    }

    /// Get bin center sizes.
    pub fn bin_sizes(&self) -> &[f64] {
        &self.bin_sizes
    }
}
