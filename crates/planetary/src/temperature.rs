//! Planetary temperature properties and calculations
//!
//! This module provides temperature calculations including:
//! - Equilibrium (blackbody) temperature from stellar radiation
//! - Effective surface temperature with greenhouse and internal heating effects
//!
//! # References
//! - Wallace & Hobbs (2006) - "Atmospheric Science: An Introductory Survey"

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::atmosphere::AtmosphereType;
use crate::planet_class::PlanetClass;

/// Planetary temperature properties
///
/// Represents both the equilibrium (blackbody) temperature from stellar radiation
/// and the effective surface temperature including atmospheric and internal effects.
///
/// # Fields
/// * `equilibrium` - Equilibrium temperature from stellar radiation balance (Kelvin)
/// * `effective` - Effective surface temperature including greenhouse effect and internal heating (Kelvin)
///
/// # Examples
/// ```
/// use planetary::temperature::Temperature;
/// use planetary::planet_class::PlanetClass;
/// use planetary::atmosphere::AtmosphereType;
///
/// // Earth: equilibrium ~255K, but greenhouse effect raises effective temp to ~288K
/// let earth_temp = Temperature::calculate(
///     1.0,
///     255.0,
///     &PlanetClass::Compact,
///     &AtmosphereType::NitrogenOxygen
/// );
/// assert_eq!(earth_temp.equilibrium, 255.0);
/// assert!((earth_temp.effective - 288.0).abs() < 5.0);
/// ```
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct Temperature {
    /// Equilibrium (blackbody) temperature from stellar radiation (Kelvin)
    pub equilibrium: f64,
    /// Effective surface temperature including greenhouse and internal heating (Kelvin)
    pub effective: f64,
}

impl Temperature {
    /// Create a new Temperature with both values
    pub fn new(equilibrium: f64, effective: f64) -> Self {
        Self {
            equilibrium,
            effective,
        }
    }

    /// Calculate effective temperature including greenhouse and internal heating effects
    ///
    /// Starts with the equilibrium (blackbody) temperature and adds contributions from:
    /// - Greenhouse effect based on atmospheric composition
    /// - Internal heating from radioactive decay and gravitational contraction
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `equilibrium` - Equilibrium (blackbody) temperature in Kelvin
    /// * `class` - Physical planet classification
    /// * `atmosphere` - Atmospheric composition
    ///
    /// # Returns
    /// Temperature struct with both equilibrium and effective temperatures
    ///
    /// # Examples
    /// ```
    /// use planetary::temperature::Temperature;
    /// use planetary::planet_class::PlanetClass;
    /// use planetary::atmosphere::AtmosphereType;
    ///
    /// // Earth: equilibrium ~255K, effective ~288K due to greenhouse effect
    /// let temp = Temperature::calculate(
    ///     1.0,
    ///     255.0,
    ///     &PlanetClass::Compact,
    ///     &AtmosphereType::NitrogenOxygen
    /// );
    /// assert!((temp.effective - 288.0).abs() < 5.0);
    /// ```
    pub fn calculate(
        mass_earth: f64,
        equilibrium: f64,
        class: &PlanetClass,
        atmosphere: &AtmosphereType,
    ) -> Self {
        // Get greenhouse effect from atmosphere
        let greenhouse_effect = atmosphere.greenhouse_effect();

        // Calculate internal heating based on planet class and mass
        let internal_heating = Self::calculate_internal_heating(mass_earth, class);

        let effective = (equilibrium + greenhouse_effect + internal_heating).round();

        Temperature {
            equilibrium,
            effective,
        }
    }

    /// Calculate internal heating contribution based on mass and class
    fn calculate_internal_heating(mass_earth: f64, class: &PlanetClass) -> f64 {
        match class {
            // Compact planets have negligible internal heating
            // Earth's internal heat contributes only ~0.1K to surface temperature
            PlanetClass::Compact => {
                if mass_earth < 0.1 {
                    0.0 // Dwarf planets/small bodies: negligible
                } else {
                    0.1 // Terrestrial planets: minimal
                }
            }

            // Transitional planets: slightly more internal heat
            // Larger cores retain more primordial heat
            PlanetClass::Transitional => 0.5 + mass_earth.powf(0.3),

            // Volatile planets (ice giants): moderate internal heating
            // Uranus: ~59K effective vs ~58K equilibrium (minimal excess)
            // Neptune: ~59K effective vs ~47K equilibrium (significant excess)
            PlanetClass::Volatile => 5.0 * (mass_earth / 15.0).powf(0.5),

            // Giant planets: significant internal heating
            // Jupiter's internal heat ~ 100K contribution
            // Scales with mass relative to Jupiter (318 M⊕)
            PlanetClass::Giant => 100.0 * (mass_earth / 318.0).powf(0.3),
        }
    }

    /// Calculate equilibrium (blackbody) temperature from stellar radiation
    ///
    /// Uses the Stefan-Boltzmann law to compute the temperature a planet would reach
    /// from stellar radiation alone, without atmospheric or internal heating effects.
    ///
    /// # Arguments
    /// * `semi_major_axis_au` - Semi-major axis in AU
    /// * `stellar_luminosity` - Stellar luminosity in solar luminosities (L☉)
    /// * `albedo` - Bond albedo (0.0-1.0), fraction of light reflected
    ///
    /// # Returns
    /// Equilibrium temperature in Kelvin
    ///
    /// # Examples
    /// ```
    /// use planetary::temperature::Temperature;
    ///
    /// // Earth at 1 AU with albedo 0.3
    /// let temp = Temperature::equilibrium_temperature(1.0, 1.0, 0.3);
    /// assert!((temp - 255.0).abs() < 5.0);
    ///
    /// // Mars at 1.52 AU with albedo 0.25
    /// let temp = Temperature::equilibrium_temperature(1.52, 1.0, 0.25);
    /// assert!((temp - 210.0).abs() < 10.0);
    /// ```
    ///
    /// # References
    /// - Wallace & Hobbs (2006) - "Atmospheric Science: An Introductory Survey"
    pub fn equilibrium_temperature(
        semi_major_axis_au: f64,
        stellar_luminosity: f64,
        albedo: f64,
    ) -> f64 {
        const SOLAR_LUMINOSITY: f64 = 3.846e26; // Watts
        const AU_IN_METERS: f64 = 1.496e11; // Meters
        const STEFAN_BOLTZMANN: f64 = 5.67e-8; // W·m⁻²·K⁻⁴

        // Convert semi-major axis from AU to meters
        let distance_m = semi_major_axis_au * AU_IN_METERS;

        // Star's total luminosity in Watts
        let star_luminosity_watts = stellar_luminosity * SOLAR_LUMINOSITY;

        // Equilibrium temperature (in Kelvin)
        // T_eq = (L * (1 - A) / (16 * π * σ * d²))^0.25
        ((star_luminosity_watts * (1.0 - albedo))
            / (16.0 * PI * STEFAN_BOLTZMANN * distance_m.powi(2)))
        .powf(0.25)
        .round()
    }

    /// Estimate albedo based on planet class and temperature
    ///
    /// Returns a reasonable default albedo for planets without specific data.
    ///
    /// # Arguments
    /// * `class` - Planet physical classification
    /// * `temperature` - Approximate temperature in Kelvin
    ///
    /// # Returns
    /// Estimated Bond albedo (0.0-1.0)
    pub fn estimate_albedo(class: &PlanetClass, temperature: f64) -> f64 {
        match class {
            PlanetClass::Giant => {
                // Gas giants: varies with cloud composition
                match temperature {
                    t if t > 1500.0 => 0.05, // Ultra-hot: dark, absorbing
                    t if t > 1000.0 => 0.10, // Hot Jupiter: low albedo
                    t if t > 500.0 => 0.30,  // Warm: intermediate
                    _ => 0.50,               // Cold: reflective ammonia clouds
                }
            }
            PlanetClass::Volatile => {
                // Ice giants: methane haze
                match temperature {
                    t if t > 500.0 => 0.20,
                    _ => 0.30, // Uranus/Neptune-like
                }
            }
            PlanetClass::Transitional => {
                // Mini-Neptunes and super-Earths: varied
                match temperature {
                    t if t > 500.0 => 0.15,
                    t if t > 300.0 => 0.30,
                    _ => 0.35,
                }
            }
            PlanetClass::Compact => {
                // Rocky planets: varies significantly
                match temperature {
                    t if t > 1000.0 => 0.10, // Lava world: dark basalt
                    t if t > 400.0 => 0.20,  // Hot desert
                    t if t > 250.0 => 0.30,  // Earth-like
                    t if t > 200.0 => 0.25,  // Mars-like
                    _ => 0.60,               // Icy: high albedo
                }
            }
        }
    }

    /// Temperature classification for this temperature
    pub fn classification(&self) -> TemperatureClass {
        TemperatureClass::classify(self.effective)
    }
}

/// Temperature classification for planets
///
/// Temperature is a modifier applied to planetary type, not a primary classification.
/// Example: "Hot Jupiter" = GasGiant with Hot temperature class
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum TemperatureClass {
    /// T < 150 K - Beyond snow line
    Cold,

    /// 150-400 K - Habitable zone range
    Temperate,

    /// 400-1000 K - Warm planets
    Warm,

    /// 1000-2000 K - Hot Jupiters, hot Neptunes
    Hot,

    /// > 2000 K - Lava worlds, ionized atmospheres
    UltraHot,
}

impl TemperatureClass {
    /// Classify temperature into categories
    ///
    /// # Arguments
    /// * `temperature` - Temperature in Kelvin
    ///
    /// # Returns
    /// Temperature classification
    pub fn classify(temperature: f64) -> Self {
        match temperature {
            t if t < 150.0 => Self::Cold,
            t if t < 400.0 => Self::Temperate,
            t if t < 1000.0 => Self::Warm,
            t if t < 2000.0 => Self::Hot,
            _ => Self::UltraHot,
        }
    }

    /// Returns the temperature range for this class in Kelvin
    pub fn range(&self) -> (f64, f64) {
        match self {
            Self::Cold => (0.0, 150.0),
            Self::Temperate => (150.0, 400.0),
            Self::Warm => (400.0, 1000.0),
            Self::Hot => (1000.0, 2000.0),
            Self::UltraHot => (2000.0, f64::INFINITY),
        }
    }

    /// Human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            Self::Cold => "Cold",
            Self::Temperate => "Temperate",
            Self::Warm => "Warm",
            Self::Hot => "Hot",
            Self::UltraHot => "Ultra-Hot",
        }
    }

    /// Whether this temperature class is potentially habitable
    pub fn potentially_habitable(&self) -> bool {
        matches!(self, Self::Temperate)
    }
}

impl std::fmt::Display for TemperatureClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}
