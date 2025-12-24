//! Atmospheric classification based on composition and pressure
//!
//! This module provides atmospheric classification for planets based on their
//! physical properties. Atmosphere type affects temperature through greenhouse
//! effects and determines surface habitability.
//!
//! # References
//! - Seager & Deming (2010) - "Exoplanet Atmospheres"
//! - Zahnle & Catling (2017) - "The Cosmic Shoreline"

use std::fmt;

use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::planet_class::PlanetClass;

/// Atmospheric classification based on composition and pressure
///
/// Classifies atmospheres by their dominant chemical species and pressure regime.
/// This affects planetary temperature through greenhouse effects and determines
/// potential habitability.
///
/// # Pressure Regimes
/// - **None/Trace**: < 10⁻¹⁰ bar (hard vacuum to exosphere)
/// - **Thin**: < 0.1 bar (Mars-like)
/// - **Medium**: 0.1-10 bar (Earth-like to Venus-like)
/// - **Thick**: > 10 bar (giant planet atmospheres)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum AtmosphereType {
    /// No atmosphere (hard vacuum)
    None,

    /// Trace atmosphere - exosphere only (<10⁻¹⁰ bar), Luna/Mercury-like
    Trace,

    // =========================================================================
    // Thin atmospheres (<0.1 bar)
    // =========================================================================
    /// Thin CO₂ atmosphere - Mars-like (~0.006 bar)
    ThinCO2,

    /// Thin nitrogen atmosphere - Pluto-like
    ThinNitrogen,

    /// Thin sulfur dioxide atmosphere - Io-like (volcanic SO₂)
    ThinSulfur,

    /// Thin oxygen atmosphere - Europa-like (radiolysis O₂)
    ThinOxygen,

    // =========================================================================
    // Medium atmospheres (0.1-10 bar)
    // =========================================================================
    /// Nitrogen-oxygen atmosphere - Earth-like N₂/O₂ dominated (~1 bar)
    NitrogenOxygen,

    /// Thick carbon dioxide atmosphere - Venus-like (~92 bar)
    ThickCO2,

    /// Water vapor dominated atmosphere - steam world
    WaterVapor,

    /// Hydrogen atmosphere - Hycean or water-rich mini-Neptune
    Hydrogen,

    // =========================================================================
    // Thick atmospheres (>10 bar)
    // =========================================================================
    /// Methane-rich atmosphere - Titan-like with N₂/CH₄
    MethaneRich,

    /// Ammonia-rich atmosphere - cold ice giant (NH₃ abundant)
    Ammonia,

    /// Hydrogen-helium atmosphere - gas giant, cool (<500K), Jupiter/Saturn-like
    HydrogenHelium,

    /// Hydrogen-methane atmosphere - ice giant, cool (<500K), Uranus/Neptune-like
    HydrogenMethane,

    /// Hydrogen-sulfide atmosphere - warm gas/ice giant (500-1000K)
    HydrogenSulfide,

    /// Hydrogen with metal oxides - hot gas giant (1000-2000K), TiO/VO clouds
    HydrogenMetalOxides,

    /// Hydrogen with ionized metals - ultra-hot gas giant (>2000K), thermal ionization
    HydrogenIonizedMetal,

    /// Unknown or unclassified atmosphere
    Unknown,
}

impl AtmosphereType {
    /// Classify atmosphere type based on planetary properties
    ///
    /// Determines atmospheric composition based on:
    /// - Atmospheric retention (escape velocity vs molecular velocities)
    /// - Planet class (gas giant, terrestrial, etc.)
    /// - Temperature (affects atmospheric chemistry and retention)
    /// - Mass (gravitational binding)
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `temperature` - Surface/equilibrium temperature in Kelvin
    /// * `class` - Planet physical class (Compact, Transitional, Volatile, Giant)
    /// * `escape_velocity` - Escape velocity in km/s
    ///
    /// # Returns
    /// Atmospheric classification enum
    ///
    /// # Examples
    /// ```
    /// use planetary::atmosphere::AtmosphereType;
    /// use planetary::planet_class::PlanetClass;
    ///
    /// // Earth-like planet: temperate, moderate mass, 11.2 km/s escape velocity
    /// let atm = AtmosphereType::classify(1.0, 288.0, &PlanetClass::Compact, 11.2);
    /// assert_eq!(atm, AtmosphereType::NitrogenOxygen);
    ///
    /// // Mars-like planet: cold, low mass, 5 km/s escape velocity
    /// let atm = AtmosphereType::classify(0.1, 210.0, &PlanetClass::Compact, 5.0);
    /// assert_eq!(atm, AtmosphereType::ThinCO2);
    /// ```
    pub fn classify(
        mass_earth: f64,
        temperature: f64,
        class: &PlanetClass,
        escape_velocity: f64,
    ) -> Self {
        // Check if body can retain an atmosphere
        // Temperature matters: cold bodies can retain atmospheres with lower escape velocity
        // (Titan has only 2.6 km/s but retains thick N2/CH4 atmosphere at 94K)
        // (Pluto has only 1.2 km/s but retains thin N2 atmosphere at 44K)
        let min_escape_velocity = match temperature {
            t if t < 50.0 => 1.0, // Ultra-cold: extremely slow molecular motion (Pluto-like)
            t if t < 150.0 => 2.0, // Very cold: slow molecular motion allows retention (Titan-like)
            t if t < 400.0 => 5.0, // Temperate: need moderate escape velocity
            _ => 7.0,             // Hot: molecules move fast, need high escape velocity
        };

        if escape_velocity < min_escape_velocity {
            // Bodies with significant mass have trace exospheres from outgassing and solar wind
            // Luna: 0.012 M⊕, ~250K, ~10^-15 bar of Ar, He, Ne
            // Mercury: 0.055 M⊕, ~440K, ~10^-14 bar of Na, K, O
            return if mass_earth > 0.001 && temperature < 500.0 {
                Self::Trace
            } else {
                Self::None
            };
        }

        match class {
            PlanetClass::Giant => Self::classify_gas_giant(temperature),
            PlanetClass::Volatile => Self::classify_ice_giant(temperature),
            PlanetClass::Transitional => {
                Self::classify_transitional(temperature, escape_velocity, mass_earth)
            }
            PlanetClass::Compact => {
                Self::classify_compact(mass_earth, temperature, escape_velocity)
            }
        }
    }

    /// Classify with additional geological context
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `temperature` - Surface/equilibrium temperature in Kelvin
    /// * `class` - Planet physical class
    /// * `escape_velocity` - Escape velocity in km/s
    /// * `is_volcanically_active` - Whether the body has active volcanism
    /// * `water_fraction` - Optional water mass fraction (for icy bodies)
    pub fn classify_with_geology(
        mass_earth: f64,
        temperature: f64,
        class: &PlanetClass,
        escape_velocity: f64,
        is_volcanically_active: bool,
        water_fraction: Option<f64>,
    ) -> Self {
        // Io-like: Small, volcanically active body with sulfur atmosphere
        // Io: 0.015 M⊕, continuous SO2 replenishment from volcanoes
        if is_volcanically_active && mass_earth < 0.05 && escape_velocity >= 2.0 {
            return Self::ThinSulfur;
        }

        // Europa-like: Small icy body with subsurface ocean
        // Europa: 0.008 M⊕, radiolysis produces tenuous O2 atmosphere
        if let Some(water) = water_fraction
            && water > 0.5
            && mass_earth < 0.02
            && temperature < 150.0
            && escape_velocity >= 1.5
        {
            return Self::ThinOxygen;
        }

        Self::classify(mass_earth, temperature, class, escape_velocity)
    }

    fn classify_gas_giant(temperature: f64) -> Self {
        match temperature {
            t if t < 500.0 => Self::HydrogenHelium,
            t if t < 1000.0 => Self::HydrogenSulfide,
            t if t < 2000.0 => Self::HydrogenMetalOxides,
            _ => Self::HydrogenIonizedMetal,
        }
    }

    fn classify_ice_giant(temperature: f64) -> Self {
        match temperature {
            t if t < 500.0 => Self::HydrogenMethane,
            t if t < 1000.0 => Self::HydrogenSulfide,
            _ => Self::HydrogenMetalOxides,
        }
    }

    fn classify_transitional(temperature: f64, escape_velocity: f64, mass_earth: f64) -> Self {
        // Transitional planets (2-5 M⊕) can be:
        // - SuperTerran (stripped cores, rocky) - thin/thick CO2
        // - MiniNeptune (retained envelope) - H/He dominated
        // - WaterWorld/Hycean - water vapor or H2

        // High water content with H2 retention → Hycean
        if escape_velocity > 10.0 && (300.0..500.0).contains(&temperature) {
            return Self::Hydrogen;
        }

        // Hot → steam world
        if temperature > 373.0 {
            return Self::WaterVapor;
        }

        // Cold with decent escape velocity → thin envelope
        if temperature < 250.0 && escape_velocity > 8.0 {
            return Self::HydrogenHelium;
        }

        // Venus-like hot super-terrestrial
        if temperature > 600.0 {
            return Self::ThickCO2;
        }

        // Default to thin CO2 for stripped cores
        if mass_earth < 3.0 {
            Self::ThinCO2
        } else {
            Self::HydrogenHelium
        }
    }

    fn classify_compact(mass_earth: f64, temperature: f64, escape_velocity: f64) -> Self {
        match (mass_earth, temperature) {
            // Pluto-like: Dwarf planet with thin nitrogen atmosphere
            (m, t) if (0.001..0.01).contains(&m) && t < 100.0 && escape_velocity >= 1.0 => {
                Self::ThinNitrogen
            }

            // Titan-like: Small moon, very cold, thick methane/nitrogen
            (m, t)
                if (0.01..=0.05).contains(&m)
                    && (80.0..120.0).contains(&t)
                    && escape_velocity >= 2.0 =>
            {
                Self::MethaneRich
            }

            // Mars-like: Small, cold, thin CO2
            (m, t) if (0.05..0.5).contains(&m) && t < 250.0 && escape_velocity < 10.0 => {
                Self::ThinCO2
            }

            // Very low mass bodies
            (m, _) if m < 0.001 => Self::None,

            // Earth-like: Temperate, nitrogen-oxygen
            (m, t) if m < 2.0 && (250.0..=350.0).contains(&t) => Self::NitrogenOxygen,

            // Venus-like: Hot, thick CO2
            (m, t) if m < 2.0 && t > 350.0 && t < 700.0 => Self::ThickCO2,

            // Hot rocky planets: Water vapor from outgassing
            (_, t) if t > 700.0 => Self::WaterVapor,

            _ => Self::None,
        }
    }

    /// Returns the greenhouse effect contribution in Kelvin
    ///
    /// This is the temperature increase from stellar equilibrium due to
    /// atmospheric heat trapping.
    pub fn greenhouse_effect(&self) -> f64 {
        match self {
            Self::None | Self::Trace => 0.0,
            Self::ThinCO2 => 5.0,
            Self::ThinNitrogen | Self::ThinSulfur | Self::ThinOxygen => 2.0,
            Self::NitrogenOxygen => 33.0, // Earth's greenhouse effect
            Self::ThickCO2 => 500.0,      // Venus-like extreme greenhouse
            Self::WaterVapor => 60.0,     // Strong H2O greenhouse
            Self::Hydrogen => 50.0,       // H2 greenhouse (can be stronger for Hycean)
            Self::MethaneRich => 20.0,    // Titan-like
            Self::Ammonia => 30.0,
            Self::HydrogenHelium
            | Self::HydrogenMethane
            | Self::HydrogenSulfide
            | Self::HydrogenMetalOxides
            | Self::HydrogenIonizedMetal => 0.0, // Giants: internal heating dominates
            Self::Unknown => 0.0,
        }
    }

    /// Returns whether this atmosphere could support liquid water
    pub fn potentially_habitable(&self) -> bool {
        matches!(
            self,
            Self::NitrogenOxygen | Self::Hydrogen | Self::WaterVapor
        )
    }

    /// Returns whether this is a primordial H/He dominated atmosphere
    pub fn is_primordial(&self) -> bool {
        matches!(
            self,
            Self::HydrogenHelium
                | Self::HydrogenMethane
                | Self::HydrogenSulfide
                | Self::HydrogenMetalOxides
                | Self::HydrogenIonizedMetal
                | Self::Hydrogen
        )
    }

    /// Human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "None",
            Self::Trace => "Trace",
            Self::ThinCO2 => "Thin CO₂",
            Self::ThinNitrogen => "Thin N₂",
            Self::ThinSulfur => "Thin SO₂",
            Self::ThinOxygen => "Thin O₂",
            Self::NitrogenOxygen => "N₂/O₂",
            Self::ThickCO2 => "Thick CO₂",
            Self::WaterVapor => "H₂O Vapor",
            Self::Hydrogen => "H₂",
            Self::MethaneRich => "N₂/CH₄",
            Self::Ammonia => "NH₃",
            Self::HydrogenHelium => "H₂/He",
            Self::HydrogenMethane => "H₂/CH₄",
            Self::HydrogenSulfide => "H₂S",
            Self::HydrogenMetalOxides => "H₂ + TiO/VO",
            Self::HydrogenIonizedMetal => "H₂ + Ionized",
            Self::Unknown => "Unknown",
        }
    }
}

impl fmt::Display for AtmosphereType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}
