//! Geological activity classification
//!
//! This module provides volcanic and tectonic activity classification
//! based on thermal state and planetary properties.
//!
//! # References
//! - Wilson & Head (1994) - "Mars review of volcanic eruption theory"
//! - Kieffer et al. (2006) - "Io's volcanic activity"
//! - Korenaga (2013) - "Initiation and evolution of plate tectonics on Earth"

use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::interior::ThermalState;
use crate::planet_class::PlanetClass;

/// Level of volcanic activity
///
/// Volcanism is driven by internal heat and provides insights into
/// planetary thermal evolution.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum VolcanismLevel {
    /// No volcanic activity (geologically dead)
    None,

    /// Ancient volcanism (>1 Gyr ago)
    Ancient,

    /// Active volcanism (present-day)
    Active,

    /// Hyperactive volcanism (tidally heated or very young)
    Hyperactive,

    /// Cryovolcanism (ice volcanism on cold bodies)
    Cryovolcanic,
}

/// Type of plate tectonics
///
/// Plate tectonics is a surface expression of mantle convection,
/// crucial for atmospheric regulation and habitability.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum TectonicRegime {
    /// No solid surface (gas/ice giants)
    None,

    /// Stagnant lid (single-plate, no motion)
    StagnantLid,

    /// Episodic lid overturn (Venus-like)
    EpisodicOverturn,

    /// Active plate tectonics (Earth-like)
    ActivePlates,

    /// Unknown/Uncertain
    Unknown,
}

/// Geological activity summary
///
/// Combines volcanism, tectonics, and surface age estimates.
///
/// # Example
/// ```
/// use planetary::geology::GeologicalActivity;
/// use planetary::interior::ThermalState;
/// use planetary::planet_class::PlanetClass;
///
/// // Earth-like geology
/// let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
/// let geology = GeologicalActivity::calculate(
///     &thermal,
///     1.0,    // mass
///     4.5,    // age
///     288.0,  // surface temp
///     &PlanetClass::Compact,
///     false,  // not water-rich
/// );
///
/// assert_eq!(geology.volcanism, planetary::geology::VolcanismLevel::Active);
/// ```
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct GeologicalActivity {
    /// Level of volcanic activity
    pub volcanism: VolcanismLevel,

    /// Type of tectonic regime
    pub tectonics: TectonicRegime,

    /// Estimated mean surface age in Gyr
    pub surface_age: f64,

    /// Estimated crustal thickness in km
    pub crustal_thickness: f64,
}

impl GeologicalActivity {
    /// Calculate geological activity from thermal state and properties
    ///
    /// # Arguments
    /// * `thermal` - Thermal state (core temp, heat flow, etc.)
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `age_gyr` - Stellar age in Gyr
    /// * `surface_temp` - Surface temperature in K
    /// * `class` - Planet physical class
    /// * `is_water_rich` - Whether planet has significant water content
    pub fn calculate(
        thermal: &ThermalState,
        mass_earth: f64,
        age_gyr: f64,
        surface_temp: f64,
        class: &PlanetClass,
        is_water_rich: bool,
    ) -> Self {
        let volcanism =
            Self::classify_volcanism(thermal, mass_earth, surface_temp, class, is_water_rich);
        let tectonics =
            Self::classify_tectonics(thermal, mass_earth, surface_temp, class, is_water_rich);
        let surface_age = Self::estimate_surface_age(thermal, mass_earth, age_gyr, &volcanism);
        let crustal_thickness = Self::estimate_crustal_thickness(mass_earth, age_gyr, class);

        Self {
            volcanism,
            tectonics,
            surface_age,
            crustal_thickness,
        }
    }

    /// Classify volcanic activity level
    fn classify_volcanism(
        thermal: &ThermalState,
        mass_earth: f64,
        surface_temp: f64,
        class: &PlanetClass,
        is_water_rich: bool,
    ) -> VolcanismLevel {
        // No volcanism for gas/ice giants (no accessible solid surface)
        if matches!(class, PlanetClass::Giant | PlanetClass::Volatile) {
            return VolcanismLevel::None;
        }

        // Cryovolcanism for cold water-rich worlds
        if is_water_rich && surface_temp < 250.0 {
            return VolcanismLevel::Cryovolcanic;
        }

        // Cold distant worlds are geologically dead unless tidally heated
        if surface_temp < 150.0 {
            if thermal.tidal_power > 10.0 * thermal.radiogenic_power && thermal.heat_flow > 0.1 {
                return VolcanismLevel::Hyperactive;
            } else {
                return VolcanismLevel::None;
            }
        }

        // Heat flow thresholds (W/m²)
        // Earth: ~0.087, Mars: ~0.016, Io: ~2.5
        match thermal.heat_flow {
            f if f > 1.0 => VolcanismLevel::Hyperactive,
            f if f > 0.04 => {
                if thermal.tidal_power > thermal.radiogenic_power {
                    VolcanismLevel::Hyperactive
                } else {
                    VolcanismLevel::Active
                }
            }
            f if f > 0.015 => {
                if mass_earth > 0.5 {
                    VolcanismLevel::Active
                } else {
                    VolcanismLevel::Ancient
                }
            }
            f if f > 0.005 => {
                if mass_earth > 0.05 {
                    VolcanismLevel::Ancient
                } else {
                    VolcanismLevel::None
                }
            }
            _ => VolcanismLevel::None,
        }
    }

    /// Classify tectonic regime
    fn classify_tectonics(
        thermal: &ThermalState,
        mass_earth: f64,
        surface_temp: f64,
        class: &PlanetClass,
        is_water_rich: bool,
    ) -> TectonicRegime {
        // No tectonics for gas/ice giants
        if matches!(class, PlanetClass::Giant | PlanetClass::Volatile) {
            return TectonicRegime::None;
        }

        // Water-rich worlds have ice shell dynamics, not traditional tectonics
        if is_water_rich && surface_temp < 273.0 {
            return TectonicRegime::StagnantLid;
        }

        // Cold distant worlds
        if surface_temp < 150.0 {
            if thermal.tidal_power > 10.0 * thermal.radiogenic_power && thermal.heat_flow > 0.1 {
                return TectonicRegime::EpisodicOverturn;
            } else {
                return TectonicRegime::StagnantLid;
            }
        }

        // Plate tectonics requires:
        // 1. High heat flow (mantle convection)
        // 2. Water (weakens lithosphere)
        // 3. Not too thick crust
        let has_high_heat_flow = thermal.heat_flow > 0.06;
        let has_liquid_water = (273.0..373.0).contains(&surface_temp);
        let right_size = mass_earth > 0.5 && mass_earth < 2.0;

        if has_high_heat_flow && has_liquid_water && right_size {
            TectonicRegime::ActivePlates
        } else if thermal.heat_flow > 0.1 && mass_earth > 0.8 {
            TectonicRegime::EpisodicOverturn
        } else {
            TectonicRegime::StagnantLid
        }
    }

    /// Estimate surface age from resurfacing rate
    fn estimate_surface_age(
        thermal: &ThermalState,
        mass_earth: f64,
        planetary_age: f64,
        volcanism: &VolcanismLevel,
    ) -> f64 {
        let resurfacing_fraction = match volcanism {
            VolcanismLevel::Hyperactive => 0.01,
            VolcanismLevel::Active => 0.1,
            VolcanismLevel::Cryovolcanic => {
                if thermal.tidal_power > 1.0e13 {
                    0.1
                } else {
                    0.8
                }
            }
            VolcanismLevel::Ancient => 0.9,
            VolcanismLevel::None => 1.0,
        };

        let mass_factor = if mass_earth < 0.5 { 1.2 } else { 1.0 };
        (planetary_age * resurfacing_fraction * mass_factor).min(planetary_age)
    }

    /// Estimate crustal thickness
    fn estimate_crustal_thickness(mass_earth: f64, age_gyr: f64, class: &PlanetClass) -> f64 {
        // Gas/ice giants don't have solid crust
        if matches!(class, PlanetClass::Giant | PlanetClass::Volatile) {
            return 0.0;
        }

        let base_thickness = match mass_earth {
            m if m < 0.1 => 20.0,
            m if m < 0.5 => 30.0 + (m - 0.1) * 50.0,
            m if m < 2.0 => 35.0 + (m - 0.5) * 10.0,
            _ => 50.0,
        };

        let age_factor = 1.0 + (age_gyr / 5.0) * 0.3;
        base_thickness * age_factor
    }

    /// Check if planet is geologically active
    pub fn is_active(&self) -> bool {
        matches!(
            self.volcanism,
            VolcanismLevel::Active | VolcanismLevel::Hyperactive | VolcanismLevel::Cryovolcanic
        )
    }

    /// Check if planet has plate tectonics
    pub fn has_plate_tectonics(&self) -> bool {
        matches!(self.tectonics, TectonicRegime::ActivePlates)
    }
}

impl VolcanismLevel {
    /// Human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "None",
            Self::Ancient => "Ancient",
            Self::Active => "Active",
            Self::Hyperactive => "Hyperactive",
            Self::Cryovolcanic => "Cryovolcanic",
        }
    }
}

impl TectonicRegime {
    /// Human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "None",
            Self::StagnantLid => "Stagnant Lid",
            Self::EpisodicOverturn => "Episodic Overturn",
            Self::ActivePlates => "Active Plates",
            Self::Unknown => "Unknown",
        }
    }
}

impl std::fmt::Display for VolcanismLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl std::fmt::Display for TectonicRegime {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}
