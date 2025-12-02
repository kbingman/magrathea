//! Size distributions for particle populations.
//!
//! Represents how particle mass is distributed across sizes.
//! Used for statistical treatment of dust, pebbles, and boulders.
//!
//! # Variants
//!
//! - **PowerLaw**: Analytical dn/ds ∝ s^(-q), efficient for initialization
//! - **Monodisperse**: Single size, useful for testing
//! - **Binned**: Discretized, required for coagulation (future)
//!
//! # Physics
//!
//! The MRN distribution (Mathis, Rumpl & Nordsieck 1977) has q = 3.5,
//! which is also the collisional equilibrium exponent (Dohnanyi 1969).
//! For q < 4, most mass is in the largest particles.
//! For q > 3, most surface area is in the smallest particles.

use std::f64::consts::PI;

use units::{Density, Length, Mass};

/// Distribution of particle sizes in a population.
///
/// All variants track the material density (internal density of particles)
/// to enable mass calculations from size.
#[derive(Debug, Clone)]
pub enum SizeDistribution {
    /// Power-law size distribution: dn/ds ∝ s^(-q)
    ///
    /// The number of particles per unit size interval follows a power law.
    /// Common values:
    /// - q = 3.5: MRN (interstellar), collisional equilibrium
    /// - q = 2.5: top-heavy, most mass in largest sizes
    PowerLaw {
        /// Minimum particle size (cm)
        s_min: f64,
        /// Maximum particle size (cm)
        s_max: f64,
        /// Power-law exponent q (dn/ds ∝ s^(-q))
        exponent: f64,
        /// Total mass in distribution (g)
        total_mass: f64,
        /// Internal material density (g/cm³)
        material_density: f64,
    },

    /// Single particle size (delta function).
    ///
    /// All particles have exactly the same size. Useful for:
    /// - Testing and validation
    /// - Simple analytical models
    /// - Representing a narrow size peak
    Monodisperse {
        /// Particle size (cm)
        size: f64,
        /// Total mass (g)
        total_mass: f64,
        /// Internal material density (g/cm³)
        material_density: f64,
    },

    /// Discretized size distribution in logarithmic bins.
    ///
    /// Required for coagulation/fragmentation calculations.
    /// Mass is tracked per bin rather than analytically.
    Binned {
        /// Bin edges in size space (cm), length N+1 for N bins
        bin_edges: Vec<f64>,
        /// Mass in each bin (g), length N
        mass_per_bin: Vec<f64>,
        /// Internal material density (g/cm³)
        material_density: f64,
    },
}

impl SizeDistribution {
    // =========================================================================
    // Constructors
    // =========================================================================

    /// Create a power-law size distribution.
    ///
    /// # Arguments
    /// * `s_min` - Minimum particle size
    /// * `s_max` - Maximum particle size
    /// * `exponent` - Power-law exponent q (dn/ds ∝ s^(-q))
    /// * `total_mass` - Total mass in distribution
    /// * `material_density` - Internal density of particles
    ///
    /// # Panics
    /// Panics if s_min >= s_max or if any value is non-positive.
    pub fn power_law(
        s_min: Length,
        s_max: Length,
        exponent: f64,
        total_mass: Mass,
        material_density: Density,
    ) -> Self {
        let s_min_cm = s_min.to_cm();
        let s_max_cm = s_max.to_cm();
        let total_mass_g = total_mass.to_grams();
        let material_density_cgs = material_density.to_grams_per_cm3();

        assert!(s_min_cm > 0.0, "s_min must be positive");
        assert!(s_max_cm > s_min_cm, "s_max must be greater than s_min");
        assert!(total_mass_g >= 0.0, "total_mass must be non-negative");
        assert!(
            material_density_cgs > 0.0,
            "material_density must be positive"
        );

        Self::PowerLaw {
            s_min: s_min_cm,
            s_max: s_max_cm,
            exponent,
            total_mass: total_mass_g,
            material_density: material_density_cgs,
        }
    }

    /// Create an MRN-like distribution (q = 3.5) from 0.1 μm to 1 μm.
    ///
    /// The MRN distribution (Mathis, Rumpl & Nordsieck 1977) describes
    /// interstellar dust grains. It's a reasonable starting point for
    /// disk dust before growth occurs.
    ///
    /// # Arguments
    /// * `total_mass` - Total mass in distribution
    /// * `material_density` - Internal density of particles
    pub fn mrn(total_mass: Mass, material_density: Density) -> Self {
        Self::power_law(
            Length::from_cm(1e-5), // 0.1 μm
            Length::from_cm(1e-4), // 1 μm
            3.5,                   // MRN exponent
            total_mass,
            material_density,
        )
    }

    /// Create a monodisperse (single-size) distribution.
    ///
    /// # Arguments
    /// * `size` - Particle size
    /// * `total_mass` - Total mass
    /// * `material_density` - Internal density of particles
    pub fn monodisperse(size: Length, total_mass: Mass, material_density: Density) -> Self {
        let size_cm = size.to_cm();
        let total_mass_g = total_mass.to_grams();
        let material_density_cgs = material_density.to_grams_per_cm3();

        assert!(size_cm > 0.0, "size must be positive");
        assert!(total_mass_g >= 0.0, "total_mass must be non-negative");
        assert!(
            material_density_cgs > 0.0,
            "material_density must be positive"
        );

        Self::Monodisperse {
            size: size_cm,
            total_mass: total_mass_g,
            material_density: material_density_cgs,
        }
    }

    /// Create a binned distribution with logarithmically-spaced bins.
    ///
    /// Initially empty (zero mass in all bins). Use `add_mass_to_bin()`
    /// or construct from another distribution with `to_binned()`.
    ///
    /// # Arguments
    /// * `s_min` - Minimum size
    /// * `s_max` - Maximum size
    /// * `n_bins` - Number of size bins
    /// * `material_density` - Internal density of particles
    pub fn binned_empty(
        s_min: Length,
        s_max: Length,
        n_bins: usize,
        material_density: Density,
    ) -> Self {
        let s_min_cm = s_min.to_cm();
        let s_max_cm = s_max.to_cm();
        let material_density_cgs = material_density.to_grams_per_cm3();

        assert!(s_min_cm > 0.0, "s_min must be positive");
        assert!(s_max_cm > s_min_cm, "s_max must be greater than s_min");
        assert!(n_bins >= 1, "need at least 1 bin");
        assert!(
            material_density_cgs > 0.0,
            "material_density must be positive"
        );

        let log_min = s_min_cm.ln();
        let log_max = s_max_cm.ln();

        let bin_edges: Vec<f64> = (0..=n_bins)
            .map(|i| {
                let frac = i as f64 / n_bins as f64;
                (log_min + frac * (log_max - log_min)).exp()
            })
            .collect();

        let mass_per_bin = vec![0.0; n_bins];

        Self::Binned {
            bin_edges,
            mass_per_bin,
            material_density: material_density_cgs,
        }
    }

    // =========================================================================
    // Queries
    // =========================================================================

    /// Total mass in the distribution.
    pub fn total_mass(&self) -> Mass {
        let mass_g = match self {
            Self::PowerLaw { total_mass, .. } => *total_mass,
            Self::Monodisperse { total_mass, .. } => *total_mass,
            Self::Binned { mass_per_bin, .. } => mass_per_bin.iter().sum(),
        };
        Mass::from_grams(mass_g)
    }

    /// Material density of particles.
    pub fn material_density(&self) -> Density {
        let rho = match self {
            Self::PowerLaw {
                material_density, ..
            } => *material_density,
            Self::Monodisperse {
                material_density, ..
            } => *material_density,
            Self::Binned {
                material_density, ..
            } => *material_density,
        };
        Density::from_grams_per_cm3(rho)
    }

    /// Minimum particle size.
    pub fn min_size(&self) -> Length {
        let s = match self {
            Self::PowerLaw { s_min, .. } => *s_min,
            Self::Monodisperse { size, .. } => *size,
            Self::Binned { bin_edges, .. } => bin_edges[0],
        };
        Length::from_cm(s)
    }

    /// Maximum particle size.
    pub fn max_size(&self) -> Length {
        let s = match self {
            Self::PowerLaw { s_max, .. } => *s_max,
            Self::Monodisperse { size, .. } => *size,
            Self::Binned { bin_edges, .. } => *bin_edges.last().unwrap(),
        };
        Length::from_cm(s)
    }

    /// Mass-weighted mean size.
    ///
    /// ⟨s⟩ = ∫ s × (dm/ds) ds / M_tot
    ///
    /// For a power-law with exponent q:
    /// - q < 4: mean is close to s_max (mass dominated by large particles)
    /// - q > 4: mean is close to s_min (mass dominated by small particles)
    pub fn mean_size(&self) -> Length {
        let s = match self {
            Self::PowerLaw {
                s_min,
                s_max,
                exponent,
                ..
            } => {
                // dm/ds ∝ s³ × s^(-q) = s^(3-q)
                // ⟨s⟩ = ∫ s × s^(3-q) ds / ∫ s^(3-q) ds
                //     = ∫ s^(4-q) ds / ∫ s^(3-q) ds
                let q = *exponent;

                if (q - 4.0).abs() < 1e-10 {
                    // q = 4: special case, dm/ds ∝ 1/s
                    // ⟨s⟩ = ∫ 1 ds / ∫ 1/s ds = (s_max - s_min) / ln(s_max/s_min)
                    (s_max - s_min) / (s_max / s_min).ln()
                } else if (q - 5.0).abs() < 1e-10 {
                    // q = 5: special case for numerator
                    let denom = (s_max.powf(4.0 - q) - s_min.powf(4.0 - q)) / (4.0 - q);
                    let numer = (s_max / s_min).ln();
                    numer / denom
                } else {
                    let numer = (s_max.powf(5.0 - q) - s_min.powf(5.0 - q)) / (5.0 - q);
                    let denom = (s_max.powf(4.0 - q) - s_min.powf(4.0 - q)) / (4.0 - q);
                    numer / denom
                }
            }
            Self::Monodisperse { size, .. } => *size,
            Self::Binned {
                bin_edges,
                mass_per_bin,
                ..
            } => {
                let total: f64 = mass_per_bin.iter().sum();
                if total == 0.0 {
                    return Length::from_cm((bin_edges[0] * bin_edges.last().unwrap()).sqrt());
                }

                let weighted_sum: f64 = bin_edges
                    .windows(2)
                    .zip(mass_per_bin.iter())
                    .map(|(edges, &mass)| {
                        let s_center = (edges[0] * edges[1]).sqrt(); // geometric mean
                        s_center * mass
                    })
                    .sum();

                weighted_sum / total
            }
        };
        Length::from_cm(s)
    }

    /// Number of particles in the distribution.
    ///
    /// N = M_tot / m_particle where m_particle = (4π/3) ρ_m s³
    ///
    /// For power-law, this requires integration over the size distribution.
    pub fn total_number(&self) -> f64 {
        match self {
            Self::PowerLaw {
                s_min,
                s_max,
                exponent,
                total_mass,
                material_density,
            } => {
                // M_tot = ∫ m(s) dn = ∫ (4π/3) ρ_m s³ × C s^(-q) ds
                //       = (4π/3) ρ_m C × ∫ s^(3-q) ds
                // N_tot = ∫ C s^(-q) ds
                //
                // So: N_tot / M_tot = [∫ s^(-q) ds] / [(4π/3) ρ_m ∫ s^(3-q) ds]
                let q = *exponent;
                let rho = *material_density;

                let mass_integral = if (q - 4.0).abs() < 1e-10 {
                    (s_max / s_min).ln()
                } else {
                    (s_max.powf(4.0 - q) - s_min.powf(4.0 - q)) / (4.0 - q)
                };

                let number_integral = if (q - 1.0).abs() < 1e-10 {
                    (s_max / s_min).ln()
                } else {
                    (s_max.powf(1.0 - q) - s_min.powf(1.0 - q)) / (1.0 - q)
                };

                let n_over_m = number_integral / ((4.0 / 3.0) * PI * rho * mass_integral);
                total_mass * n_over_m
            }
            Self::Monodisperse {
                size,
                total_mass,
                material_density,
            } => {
                let m_particle = (4.0 / 3.0) * PI * material_density * size.powi(3);
                total_mass / m_particle
            }
            Self::Binned {
                bin_edges,
                mass_per_bin,
                material_density,
            } => bin_edges
                .windows(2)
                .zip(mass_per_bin.iter())
                .map(|(edges, &mass)| {
                    let s_center = (edges[0] * edges[1]).sqrt();
                    let m_particle = (4.0 / 3.0) * PI * material_density * s_center.powi(3);
                    mass / m_particle
                })
                .sum(),
        }
    }

    /// Number of size bins (for Binned variant).
    ///
    /// Returns 1 for Monodisperse, 0 for PowerLaw (must convert first).
    pub fn n_bins(&self) -> usize {
        match self {
            Self::PowerLaw { .. } => 0,
            Self::Monodisperse { .. } => 1,
            Self::Binned { mass_per_bin, .. } => mass_per_bin.len(),
        }
    }

    /// Iterator over (bin_center, mass) pairs for binned distributions.
    ///
    /// For Binned variant, yields the geometric center of each bin and its mass.
    /// For Monodisperse, yields the single (size, mass) pair.
    /// For PowerLaw, returns empty iterator (convert to binned first).
    ///
    /// # Example
    /// ```
    /// use stellar_forge::SizeDistribution;
    /// use units::{Density, Length, Mass};
    ///
    /// let dist = SizeDistribution::mrn(
    ///     Mass::from_grams(1000.0),
    ///     Density::from_grams_per_cm3(3.0),
    /// ).to_binned(10);
    ///
    /// for (size, mass) in dist.bins() {
    ///     println!("size = {} cm, mass = {} g", size.to_cm(), mass.to_grams());
    /// }
    /// ```
    pub fn bins(&self) -> Vec<(Length, Mass)> {
        match self {
            Self::PowerLaw { .. } => vec![],
            Self::Monodisperse {
                size, total_mass, ..
            } => {
                vec![(Length::from_cm(*size), Mass::from_grams(*total_mass))]
            }
            Self::Binned {
                bin_edges,
                mass_per_bin,
                ..
            } => bin_edges
                .windows(2)
                .zip(mass_per_bin.iter())
                .map(|(edges, &mass)| {
                    let center = (edges[0] * edges[1]).sqrt();
                    (Length::from_cm(center), Mass::from_grams(mass))
                })
                .collect(),
        }
    }

    /// Convert to a binned representation with the given number of bins.
    ///
    /// Preserves total mass. Useful for coagulation calculations.
    pub fn to_binned(&self, n_bins: usize) -> Self {
        let s_min = self.min_size().to_cm();
        let s_max = self.max_size().to_cm();
        let material_density = self.material_density().to_grams_per_cm3();

        match self {
            Self::Binned { .. } => self.clone(),
            Self::Monodisperse {
                size, total_mass, ..
            } => {
                // Put all mass in the bin containing this size
                let mut result = Self::binned_empty(
                    Length::from_cm(s_min * 0.5), // extend range slightly
                    Length::from_cm(s_max * 2.0),
                    n_bins,
                    Density::from_grams_per_cm3(material_density),
                );

                if let Self::Binned {
                    bin_edges,
                    mass_per_bin,
                    ..
                } = &mut result
                {
                    // Find the bin containing this size
                    let bin_idx = bin_edges
                        .windows(2)
                        .position(|e| *size >= e[0] && *size < e[1])
                        .unwrap_or(n_bins - 1);
                    mass_per_bin[bin_idx] = *total_mass;
                }

                result
            }
            Self::PowerLaw {
                s_min,
                s_max,
                exponent,
                total_mass,
                material_density,
            } => {
                let log_min = s_min.ln();
                let log_max = s_max.ln();
                let q = *exponent;

                let bin_edges: Vec<f64> = (0..=n_bins)
                    .map(|i| {
                        let frac = i as f64 / n_bins as f64;
                        (log_min + frac * (log_max - log_min)).exp()
                    })
                    .collect();

                // Mass in each bin: M_bin = ∫_{s1}^{s2} (dm/ds) ds
                // where dm/ds ∝ s^(3-q)
                let mass_integral_total = if (q - 4.0).abs() < 1e-10 {
                    (s_max / s_min).ln()
                } else {
                    (s_max.powf(4.0 - q) - s_min.powf(4.0 - q)) / (4.0 - q)
                };

                let mass_per_bin: Vec<f64> = bin_edges
                    .windows(2)
                    .map(|edges| {
                        let s1 = edges[0];
                        let s2 = edges[1];

                        let mass_integral_bin = if (q - 4.0).abs() < 1e-10 {
                            (s2 / s1).ln()
                        } else {
                            (s2.powf(4.0 - q) - s1.powf(4.0 - q)) / (4.0 - q)
                        };

                        total_mass * mass_integral_bin / mass_integral_total
                    })
                    .collect();

                Self::Binned {
                    bin_edges,
                    mass_per_bin,
                    material_density: *material_density,
                }
            }
        }
    }
}
