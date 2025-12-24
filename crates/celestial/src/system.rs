//! Planetary system structure

use serde::{Deserialize, Serialize};
use stellar::StellarObject;
use units::{EARTH_MASS_G, SOLAR_MASS_G};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::architecture::SystemArchitecture;
use crate::metadata::SystemMetadata;
use planetary::planet::Planet;
use planetary::planet_class::PlanetClass;

/// Earth masses per solar mass (M☉/M⊕)
const EARTH_MASSES_PER_SOLAR: f64 = SOLAR_MASS_G / EARTH_MASS_G;

/// A complete planetary system with one or more stellar hosts
///
/// This is the unified output format for all system generation approaches:
/// statistical sampling, stellar-forge formation simulation, and manual construction.
///
/// # Examples
///
/// ```
/// use celestial::{PlanetarySystem, SystemArchitecture, SystemMetadata, GenerationMethod};
/// use stellar::StellarObject;
///
/// // Systems are typically created via generation functions,
/// // but can also be constructed directly for testing or manual input.
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct PlanetarySystem {
    /// Stellar host(s) - at least one required
    ///
    /// For single-star systems, this contains one element.
    /// For binary systems, contains two elements (primary first).
    /// For hierarchical systems, elements are ordered by hierarchy.
    pub stars: Vec<StellarObject>,

    /// Planets in the system, sorted by semi-major axis
    pub planets: Vec<Planet>,

    /// System metadata (generation info, architecture, identification)
    pub metadata: SystemMetadata,
}

impl PlanetarySystem {
    /// Create a new planetary system
    ///
    /// Planets are automatically sorted by semi-major axis.
    ///
    /// # Panics
    /// Panics if `stars` is empty - every system must have at least one star.
    pub fn new(
        stars: Vec<StellarObject>,
        mut planets: Vec<Planet>,
        metadata: SystemMetadata,
    ) -> Self {
        assert!(
            !stars.is_empty(),
            "PlanetarySystem must have at least one star"
        );

        planets.sort_by(|a, b| {
            a.semi_major_axis
                .to_au()
                .partial_cmp(&b.semi_major_axis.to_au())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        Self {
            stars,
            planets,
            metadata,
        }
    }

    /// Returns the primary (first) star in the system
    ///
    /// For single-star systems, this is the only star.
    /// For multi-star systems, this is typically the most massive component.
    pub fn primary_star(&self) -> &StellarObject {
        &self.stars[0]
    }

    /// Total luminosity of all stellar components (L☉)
    ///
    /// Used for habitable zone calculations in multi-star systems.
    pub fn total_luminosity(&self) -> f64 {
        self.stars.iter().map(|s| s.luminosity()).sum()
    }

    /// Effective stellar mass for orbital dynamics (M☉)
    ///
    /// For single stars, returns the star's mass.
    /// For close binaries, returns combined mass.
    pub fn effective_mass(&self) -> f64 {
        self.stars.iter().map(|s| s.mass().to_solar_masses()).sum()
    }

    /// Primary star's metallicity [Fe/H]
    pub fn metallicity(&self) -> f64 {
        self.primary_star().metallicity()
    }

    /// Primary star's spectral type as string (e.g., "G2")
    pub fn spectral_type(&self) -> String {
        self.primary_star().spectral_type_string()
    }

    /// Whether this is a multi-star system
    pub fn is_binary(&self) -> bool {
        self.stars.len() > 1
    }

    /// System architecture classification
    pub fn architecture(&self) -> SystemArchitecture {
        self.metadata.architecture
    }

    /// Check Hill stability between adjacent planet pairs
    ///
    /// Returns true if all adjacent planet pairs have sufficient separation
    /// (> 8 mutual Hill radii) for long-term stability.
    pub fn is_stable(&self) -> bool {
        if self.planets.len() < 2 {
            return true;
        }

        let stellar_mass = self.effective_mass();

        for window in self.planets.windows(2) {
            let inner = &window[0];
            let outer = &window[1];

            let m1 = inner.mass.to_earth_masses() / EARTH_MASSES_PER_SOLAR;
            let m2 = outer.mass.to_earth_masses() / EARTH_MASSES_PER_SOLAR;
            let a1 = inner.semi_major_axis.to_au();
            let a2 = outer.semi_major_axis.to_au();

            let mutual_hill =
                ((m1 + m2) / (3.0 * stellar_mass)).powf(1.0 / 3.0) * ((a1 + a2) / 2.0);
            let separation = a2 - a1;

            if separation / mutual_hill < 8.0 {
                return false;
            }
        }

        true
    }

    /// Filter planets by class
    pub fn planets_of_class(&self, class: PlanetClass) -> Vec<&Planet> {
        self.planets.iter().filter(|p| p.class == class).collect()
    }

    /// Find planets within the habitable zone
    pub fn habitable_zone_planets(&self) -> Vec<&Planet> {
        let luminosity = self.total_luminosity();
        self.planets
            .iter()
            .filter(|p| p.in_habitable_zone(luminosity))
            .collect()
    }

    /// Habitable zone boundaries for this system
    pub fn habitable_zone(&self) -> HabitableZone {
        HabitableZone::from_luminosity(self.total_luminosity())
    }

    /// Snow line distance in AU
    pub fn snow_line(&self) -> f64 {
        snow_line(self.total_luminosity())
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
