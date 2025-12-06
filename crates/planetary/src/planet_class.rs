//! Primary planet classification by physical mass regime
//!
//! Based on Wolfgang, Rogers, & Ford (2016) with refined boundaries for envelope physics.
//! This is the first tier of the two-tier classification system.

use rand::Rng;
use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

/// Classification of planet by physical mass regime
///
/// Based on Wolfgang, Rogers, & Ford (2016) with refined boundaries for envelope physics.
/// Mass-radius relationships follow distinct physical regimes determined by internal structure.
///
/// | Class        | Mass Range        | Physical Regime           | Radius Behavior              |
/// |--------------|-------------------|---------------------------|------------------------------|
/// | Rocky        | < 2 M⊕            | Self-compression          | R ∝ M^0.29, ρ increases      |
/// | Transitional | 2-5 M⊕            | Thin envelope regime      | R ∝ M^0.40, variable ρ       |
/// | Volatile     | 5-160 M⊕          | Thick envelope regime     | R ∝ M^0.56, ρ decreases      |
/// | Giant        | > 160 M⊕ (~0.5 Mj)| Electron degeneracy       | R nearly constant, ρ ∝ M     |
///
/// **Key boundaries:**
/// - **2 M⊕**: Fulton gap / radius valley - minimum mass to retain H/He envelopes
/// - **5 M⊕**: Upper end of Fulton gap - transition to sub-Neptune regime
/// - **160 M⊕**: Electron degeneracy onset - Jupiter-like behavior begins
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum PlanetClass {
    /// Rocky planets (< 2 M⊕)
    /// Self-compression dominated, cannot retain H/He envelopes
    /// Examples: Earth, Mars, Venus, rocky exoplanets
    Rocky,

    /// Transitional planets (2-5 M⊕)
    /// Thin H/He envelopes (1-10% by mass), super-Earths and mini-Neptunes
    /// Envelope retention depends on photoevaporation and formation history
    /// Examples: Some Kepler super-Earths, mini-Neptunes in Fulton gap
    Transitional,

    /// Volatile planets (5-160 M⊕)
    /// Thick H/He envelopes (10-80% by mass), sub-Neptunes and ice giants
    /// Envelope dominates structure and radius
    /// Examples: Neptune, Uranus, sub-Neptunes
    Volatile,

    /// Giant planets (> 160 M⊕, ~0.5 Jupiter masses)
    /// Electron degeneracy pressure dominates, radius nearly constant with mass
    /// Examples: Jupiter, Saturn, hot Jupiters
    Giant,
}

impl PlanetClass {
    /// Mass threshold between Rocky and Transitional regimes (2 Earth masses)
    pub const ROCKY_TRANSITIONAL_THRESHOLD: f64 = 2.0;

    /// Mass threshold between Transitional and Volatile regimes (5 Earth masses)
    pub const TRANSITIONAL_VOLATILE_THRESHOLD: f64 = 5.0;

    /// Mass threshold between Volatile and Giant regimes (~0.5 Jupiter masses)
    pub const VOLATILE_GIANT_THRESHOLD: f64 = 160.0;

    /// Classify a planet by its mass in Earth masses
    pub fn from_earth_masses(mass_earth: f64) -> Self {
        match mass_earth {
            m if m < Self::ROCKY_TRANSITIONAL_THRESHOLD => Self::Rocky,
            m if m < Self::TRANSITIONAL_VOLATILE_THRESHOLD => Self::Transitional,
            m if m < Self::VOLATILE_GIANT_THRESHOLD => Self::Volatile,
            _ => Self::Giant,
        }
    }

    /// Human-readable name for the planet class
    pub fn name(&self) -> &'static str {
        match self {
            Self::Rocky => "Rocky",
            Self::Transitional => "Transitional",
            Self::Volatile => "Volatile",
            Self::Giant => "Giant",
        }
    }

    /// Get mass-radius power law parameters for this regime
    ///
    /// Returns (coefficient, exponent, scatter_sigma) where R = coeff * M^exp (Earth units)
    ///
    /// Scatter sigma is in log10 space - a value of 0.05 means ~12% variation (10^0.05 ≈ 1.12)
    ///
    /// # References
    /// - Chen & Kipping (2017) for rocky/transitional planets
    /// - Wolfgang et al. (2016) for volatile planets
    /// - Thorngren et al. (2016) for giant planets
    pub fn mass_radius_params(&self) -> (f64, f64, f64) {
        match self {
            // Rocky: R = M^0.27 for M < 2 M⊕ (Chen & Kipping 2017)
            // Scatter ~8% (σ_log ≈ 0.035)
            Self::Rocky => (1.0, 0.27, 0.035),
            // Transitional: steeper slope as envelope physics matter
            // Scatter ~15% due to varied envelope fractions
            Self::Transitional => (1.0, 0.35, 0.06),
            // Volatile: R ∝ M^0.55 for sub-Neptunes/ice giants
            // Higher scatter ~20% due to envelope inflation variation
            Self::Volatile => (1.0, 0.55, 0.08),
            // Giant: nearly constant radius (electron degeneracy)
            // R ≈ 11 R⊕ with slight decrease at higher masses
            // Low scatter ~10% for mature giants
            Self::Giant => (11.2, 0.01, 0.04),
        }
    }

    /// Calculate radius from mass using empirical mass-radius relation
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `add_scatter` - Whether to add observational scatter
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// Radius in Earth radii
    pub fn radius_from_mass(&self, mass_earth: f64, add_scatter: bool, rng: &mut impl Rng) -> f64 {
        let (coeff, exp, sigma) = self.mass_radius_params();
        let base_radius = coeff * mass_earth.powf(exp);

        if add_scatter {
            // Scatter in log10 space - sigma of 0.05 gives ~12% variation
            let log_scatter = rng.random_range(-sigma..sigma);
            base_radius * 10_f64.powf(log_scatter)
        } else {
            base_radius
        }
    }

    /// Returns the mass range for this class in Earth masses
    pub fn mass_range(&self) -> (f64, f64) {
        match self {
            Self::Rocky => (0.0, Self::ROCKY_TRANSITIONAL_THRESHOLD),
            Self::Transitional => (
                Self::ROCKY_TRANSITIONAL_THRESHOLD,
                Self::TRANSITIONAL_VOLATILE_THRESHOLD,
            ),
            Self::Volatile => (
                Self::TRANSITIONAL_VOLATILE_THRESHOLD,
                Self::VOLATILE_GIANT_THRESHOLD,
            ),
            Self::Giant => (Self::VOLATILE_GIANT_THRESHOLD, f64::INFINITY),
        }
    }

    /// Returns whether this class can retain a primordial H/He envelope
    pub fn can_retain_envelope(&self) -> bool {
        matches!(self, Self::Transitional | Self::Volatile | Self::Giant)
    }
}

impl std::fmt::Display for PlanetClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}
