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

impl SystemArchitecture {
    pub fn sample(rng: &mut ChaChaRng, spectral_type: &str, metallicity: f64) -> Self {
        let giant_prob = 0.10 * 10.0_f64.powf(2.0 * metallicity);

        match spectral_type {
            "M" => {
                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.45 => Self::CompactMulti,
                    x if x < 0.70 => Self::Sparse,
                    x if x < 0.70 + giant_prob * 0.3 => Self::GiantDominated,
                    _ => Self::Mixed,
                }
            }
            "K" | "G" => {
                let roll: f64 = rng.random();
                match roll {
                    x if x < 0.25 => Self::CompactMulti,
                    x if x < 0.50 => Self::Mixed,
                    x if x < 0.50 + giant_prob => Self::GiantDominated,
                    _ => Self::Sparse,
                }
            }
            "F" | "A" | "B" => {
                let roll: f64 = rng.random();
                match roll {
                    x if x < giant_prob * 1.5 => Self::GiantDominated,
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
