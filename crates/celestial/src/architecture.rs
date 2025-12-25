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
    /// Architecture probabilities reflect that MOST stars have planets:
    /// - Sparse systems are RARE (< 10%), not the default
    /// - CompactMulti (TRAPPIST-1 style) is common for M/K dwarfs
    /// - Mixed (Solar System style) is common for G/K/F stars
    /// - GiantDominated scales with metallicity and stellar mass
    ///
    /// These rates are intentionally higher than Kepler detections to reflect
    /// true populations. The Solar System is probably typical, not special.
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

        // Stellar mass scaling for giant formation
        let mass_scaling = if stellar_mass < MIN_STELLAR_MASS_FOR_GIANTS {
            0.0 // No giant-dominated architecture for very low mass stars
        } else if stellar_mass < 0.5 {
            stellar_mass.powf(2.0) // M-dwarfs: steep suppression
        } else if stellar_mass < 0.9 {
            stellar_mass.powf(1.0) // K-dwarfs: linear scaling
        } else if stellar_mass < 1.2 {
            let t = (stellar_mass - 0.9) / 0.3;
            0.9 + t * 0.1 // G-dwarfs: smooth transition
        } else {
            stellar_mass.powf(1.5) // F/A stars: enhanced occurrence
        };

        let giant_factor = metallicity_boost * mass_scaling;
        let spectral_class = spectral_type.chars().next().unwrap_or('G');

        match spectral_class {
            'M' => {
                // M-dwarfs: Compact multi-planet systems dominate (TRAPPIST-1 style)
                // Sparse is RARE - most M-dwarfs have planets
                let giant_prob = (0.08 * giant_factor).min(0.15);

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.55 => Self::CompactMulti, // Dominant architecture
                    x if x < 0.85 => Self::Mixed,        // Many also have outer planets
                    x if x < 0.85 + giant_prob => Self::GiantDominated,
                    x if x < 0.95 => Self::CompactMulti, // More compact systems
                    _ => Self::Sparse,                   // Rare: ~5%
                }
            }
            'K' => {
                // K-dwarfs: Mix of compact and mixed systems
                let giant_prob = (0.10 * giant_factor).min(0.20);

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.35 => Self::CompactMulti,
                    x if x < 0.80 => Self::Mixed, // Solar System style common
                    x if x < 0.80 + giant_prob => Self::GiantDominated,
                    x if x < 0.95 => Self::Mixed, // More mixed systems
                    _ => Self::Sparse,            // Rare: ~5%
                }
            }
            'G' => {
                // G-dwarfs: Solar System style is typical
                let giant_prob = (0.08 * giant_factor).min(0.20);

                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.30 => Self::CompactMulti,
                    x if x < 0.75 => Self::Mixed, // Solar System is typical!
                    x if x < 0.75 + giant_prob => Self::GiantDominated,
                    x if x < 0.93 => Self::Mixed, // Even more mixed
                    _ => Self::Sparse,            // Rare: ~7%
                }
            }
            'F' | 'A' | 'B' => {
                // Massive stars: enhanced giant probability, mixed architectures
                let giant_prob = (0.20 * giant_factor).min(0.35);

                let roll: f64 = rng.random();
                match roll {
                    x if x < giant_prob => Self::GiantDominated,
                    x if x < 0.70 => Self::Mixed,
                    x if x < 0.85 => Self::CompactMulti,
                    _ => Self::Sparse, // Rare: ~15% (shorter disk lifetimes)
                }
            }
            _ => Self::Mixed, // Default to having planets, not sparse
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
