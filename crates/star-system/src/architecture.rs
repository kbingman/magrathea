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

        // Stellar mass scaling: P(giant) ∝ M^2.0
        // Hard cutoff below MIN_STELLAR_MASS_FOR_GIANTS
        // For 0.3 M☉: 0.09, for 0.5 M☉: 0.25, for 1.0 M☉: 1.0
        let mass_scaling = if stellar_mass < MIN_STELLAR_MASS_FOR_GIANTS {
            0.0 // No giant-dominated architecture for very low mass stars
        } else {
            stellar_mass.powf(2.0)
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
            'K' | 'G' => {
                // K/G stars have moderate base giant probability
                // At solar metallicity and mass: 0.06 * 1.0 = 6% GiantDominated
                // Combined with outer system cold giants, gives ~12-15% total
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
                // Massive stars: enhanced giant probability
                // ~8% base at solar metallicity and 1.3 M☉, scaled by mass
                // Target: 15-20% total giant occurrence for F/A stars
                let giant_prob = 0.08 * giant_factor;

                let roll: f64 = rng.random();
                match roll {
                    x if x < giant_prob => Self::GiantDominated,
                    x if x < 0.20 => Self::Mixed,
                    _ => Self::Sparse,
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
