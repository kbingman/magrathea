//! Planetary system structure and architecture

use rand::Rng;
use rand_chacha::ChaChaRng;
use serde::{Deserialize, Serialize};

use crate::planet::Planet;
use crate::planet_class::PlanetClass;

/// A complete planetary system
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PlanetarySystem {
    pub stellar_mass: f64,
    pub stellar_luminosity: f64,
    pub stellar_temperature: f64,
    pub stellar_metallicity: f64,
    pub spectral_type: String,
    pub planets: Vec<Planet>,
    pub architecture: SystemArchitecture,
}

impl PlanetarySystem {
    pub fn new(
        stellar_mass: f64,
        stellar_luminosity: f64,
        stellar_temperature: f64,
        stellar_metallicity: f64,
        spectral_type: String,
        mut planets: Vec<Planet>,
        architecture: SystemArchitecture,
    ) -> Self {
        planets.sort_by(|a, b| {
            a.semi_major_axis
                .to_au()
                .partial_cmp(&b.semi_major_axis.to_au())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        Self {
            stellar_mass,
            stellar_luminosity,
            stellar_temperature,
            stellar_metallicity,
            spectral_type,
            planets,
            architecture,
        }
    }

    /// Check Hill stability
    pub fn is_stable(&self) -> bool {
        if self.planets.len() < 2 {
            return true;
        }

        for window in self.planets.windows(2) {
            let inner = &window[0];
            let outer = &window[1];

            let m1 = inner.mass.to_earth_masses() / 332946.0;
            let m2 = outer.mass.to_earth_masses() / 332946.0;
            let a1 = inner.semi_major_axis.to_au();
            let a2 = outer.semi_major_axis.to_au();

            let mutual_hill = ((m1 + m2) / (3.0 * self.stellar_mass)).powf(1.0 / 3.0) * ((a1 + a2) / 2.0);
            let separation = a2 - a1;

            if separation / mutual_hill < 8.0 {
                return false;
            }
        }

        true
    }

    pub fn planets_of_class(&self, class: PlanetClass) -> Vec<&Planet> {
        self.planets.iter().filter(|p| p.class == class).collect()
    }

    pub fn habitable_zone_planets(&self) -> Vec<&Planet> {
        self.planets
            .iter()
            .filter(|p| p.in_habitable_zone(self.stellar_luminosity))
            .collect()
    }
}

/// System architecture classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
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

/// Snow line location in AU
pub fn snow_line(stellar_luminosity: f64) -> f64 {
    2.7 * stellar_luminosity.sqrt()
}

/// Habitable zone boundaries
pub struct HabitableZone {
    pub inner_edge: f64,
    pub outer_edge: f64,
}

impl HabitableZone {
    pub fn from_luminosity(luminosity: f64) -> Self {
        Self {
            inner_edge: (luminosity / 1.1).sqrt(),
            outer_edge: (luminosity / 0.36).sqrt(),
        }
    }
}
