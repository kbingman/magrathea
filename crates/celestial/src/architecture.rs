//! System architecture classification

use rand::Rng;
use rand_chacha::ChaChaRng;
use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

/// System architecture classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum SystemArchitecture {
    CompactMulti,
    Mixed,
    GiantDominated,
    Sparse,
}

/// Minimum stellar mass for giant-dominated architecture (solar masses)
/// Below this mass, disk mass and lifetime are insufficient for gas giant formation
/// via core accretion. Based on Laughlin+ 2004, Ida & Lin 2005.
const MIN_STELLAR_MASS_FOR_GIANTS: f64 = 0.25;

impl SystemArchitecture {
    /// Sample a system architecture based on stellar properties
    ///
    /// # Arguments
    /// * `spectral_type` - Full spectral type string (e.g., "M6", "K4", "G2")
    /// * `stellar_mass` - Stellar mass in solar masses (M☉)
    /// * `metallicity` - Stellar metallicity [Fe/H] in dex
    ///
    /// Architecture probabilities are calibrated against Kepler occurrence rates:
    /// - M-dwarfs: High compact multi-planet rate (~45%), moderate sparse
    /// - K/G-dwarfs: Balanced distribution, metallicity-dependent giant rate
    /// - F/A/B stars: Higher sparse rate, giant-dominated when metal-rich
    ///
    /// Giant-dominated architecture probability scales with stellar mass (M^2.0)
    /// to reflect the core accretion timescale problem for low-mass stars.
    /// Stars below 0.25 M☉ cannot form gas giants via core accretion.
    ///
    /// # References
    /// - Johnson et al. (2010) - Giant planet occurrence vs stellar mass
    /// - Mulders et al. (2015) - Planet occurrence around M-dwarfs
    /// - Laughlin et al. (2004) - Giant planet formation around low-mass stars
    pub fn sample(
        rng: &mut ChaChaRng,
        spectral_type: &str,
        stellar_mass: f64,
        metallicity: f64,
    ) -> Self {
        // Metallicity boost for giant formation: P(giant) ∝ 10^(2×[Fe/H])
        let metallicity_boost = 10.0_f64.powf(2.0 * metallicity);

        // Stellar mass scaling: different regimes for different mass ranges
        // - M-dwarfs: M^2.0 (steep suppression)
        // - K-dwarfs: M^1.0 (linear scaling)
        // - G-dwarfs: smooth transition to baseline
        // - F/A stars: M^1.5 (enhanced occurrence)
        let mass_scaling = if stellar_mass < MIN_STELLAR_MASS_FOR_GIANTS {
            0.0 // No giant-dominated architecture for very low mass stars
        } else if stellar_mass < 0.5 {
            // M-dwarfs: steep suppression
            stellar_mass.powf(2.0)
        } else if stellar_mass < 0.9 {
            // K-dwarfs: linear scaling
            stellar_mass.powf(1.0)
        } else if stellar_mass < 1.2 {
            // G-dwarfs: smooth transition
            let low = 0.9;
            let high = 1.0;
            let t = (stellar_mass - 0.9) / 0.3;
            low + t * (high - low)
        } else {
            // F/A stars: enhanced occurrence
            stellar_mass.powf(1.5)
        };

        // Combined giant probability factor
        let giant_factor = metallicity_boost * mass_scaling;

        // Extract the primary spectral class (first character)
        let spectral_class = spectral_type.chars().next().unwrap_or('G');

        match spectral_class {
            'M' => {
                // Base giant prob for M-dwarfs is low, scaled by mass and metallicity
                // Early M (~0.4 M☉): ~3% giant-dominated
                // Late M (~0.1 M☉): ~0.3% giant-dominated
                let giant_prob = 0.10 * giant_factor;

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.45 => Self::CompactMulti,
                    x if x < 0.70 => Self::Sparse,
                    x if x < 0.70 + giant_prob => Self::GiantDominated,
                    _ => Self::Mixed,
                }
            }
            'K' => {
                // K-dwarfs: moderate giant probability, favor Mixed architecture
                // At solar metallicity and 0.7 M☉: 0.08 * 1.0 * 0.7 = 5.6% GiantDominated
                // Combined with 40% Mixed (1.5× modifier), gives ~10-12% total
                let giant_prob = 0.08 * giant_factor;

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.20 => Self::CompactMulti,
                    x if x < 0.60 => Self::Mixed, // K-stars favor Mixed like F/A
                    x if x < 0.60 + giant_prob => Self::GiantDominated,
                    _ => Self::Sparse,
                }
            }
            'G' => {
                // G-dwarfs: balanced distribution
                // At solar metallicity and mass: 0.06 * 1.0 = 6% GiantDominated
                // Combined with 25% Mixed (1.5× modifier), gives ~10-15% total
                let giant_prob = 0.06 * giant_factor;

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.25 => Self::CompactMulti,
                    x if x < 0.50 => Self::Mixed,
                    x if x < 0.50 + giant_prob => Self::GiantDominated,
                    _ => Self::Sparse,
                }
            }
            'F' | 'A' | 'B' => {
                // Massive stars: enhanced giant probability and favor Mixed architecture
                // F/A stars have massive protoplanetary disks and long formation timescales
                // Target: 30-50% total giant occurrence for F/A stars
                // ~15% from GiantDominated + ~35% from Mixed architecture
                let giant_prob = 0.15 * giant_factor;

                let roll: f64 = rng.random();
                match roll {
                    x if x < giant_prob => Self::GiantDominated,
                    x if x < 0.50 => Self::Mixed, // F/A stars favor Solar System-like architectures
                    x if x < 0.70 => Self::Sparse,
                    _ => Self::CompactMulti, // Rare but possible
                }
            }
            _ => Self::Sparse,
        }
    }

    pub fn expected_planet_count(&self) -> (usize, usize) {
        match self {
            Self::CompactMulti => (4, 8),
            Self::Mixed => (2, 5),
            Self::GiantDominated => (1, 3),
            Self::Sparse => (0, 1),
        }
    }
}

impl std::fmt::Display for SystemArchitecture {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CompactMulti => write!(f, "Compact Multi-planet"),
            Self::Mixed => write!(f, "Mixed"),
            Self::GiantDominated => write!(f, "Giant-dominated"),
            Self::Sparse => write!(f, "Sparse"),
        }
    }
}
