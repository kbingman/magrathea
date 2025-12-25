//! Planetary interior thermal state and differentiation
//!
//! This module provides thermal state calculations for planetary bodies,
//! including internal heating sources and heat flow estimates.
//!
//! # References
//! - Schubert et al. (2001) - "Mantle Convection in the Earth and Planets"
//! - Peale et al. (1979) - "Melting of Io by tidal dissipation"

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

/// Differentiation state of planetary interior
///
/// Represents how separated a planet's internal structure is, from primitive
/// homogeneous composition to fully differentiated core-mantle-crust structure.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum DifferentiationState {
    /// Primitive, homogeneous composition (age < 10 Myr or mass < 0.1 M⊕)
    Undifferentiated,

    /// Core formation in progress (10-100 Myr)
    PartiallyDifferentiated,

    /// Fully separated into core-mantle-crust structure (age > 100 Myr)
    FullyDifferentiated,
}

/// Sources of internal heat for planetary bodies
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum HeatSource {
    /// Accretion heat (early formation only, < 100 Myr)
    Accretion,

    /// Radiogenic decay (U-238, Th-232, K-40)
    Radiogenic,

    /// Tidal heating from orbital eccentricity
    Tidal,

    /// Primordial heat (leftover from formation)
    Primordial,
}

/// Thermal state of planetary interior
///
/// Describes the current thermal energy budget and heat flow.
///
/// # Example
/// ```
/// use planetary::interior::ThermalState;
///
/// // Earth-like thermal state
/// let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
///
/// assert!(thermal.core_temperature > 4000.0); // Earth's core ~5200 K
/// assert!(thermal.heat_flow > 0.05); // Active heat loss
/// ```
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct ThermalState {
    /// Core temperature in Kelvin
    pub core_temperature: f64,

    /// Mantle temperature in Kelvin
    pub mantle_temperature: f64,

    /// Surface heat flow in W/m²
    pub heat_flow: f64,

    /// Primary heat sources (ordered by contribution)
    pub heat_sources: Vec<HeatSource>,

    /// Radiogenic heating power (W)
    pub radiogenic_power: f64,

    /// Tidal heating power (W)
    pub tidal_power: f64,
}

impl ThermalState {
    /// Calculate thermal state from planetary properties
    ///
    /// # Arguments
    /// * `mass_earth` - Planet mass in Earth masses (M⊕)
    /// * `age_gyr` - Stellar age in Gyr (proxy for planet age)
    /// * `tidal_heating_factor` - Normalized tidal heating (1.0 = Earth-Moon baseline)
    ///
    /// # Returns
    /// Thermal state with temperatures and heat sources
    ///
    /// # Example
    /// ```
    /// use planetary::interior::ThermalState;
    ///
    /// // Earth-like planet (minimal tidal heating)
    /// let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
    /// assert!(thermal.core_temperature > 4000.0);
    ///
    /// // Io-like tidally heated moon (factor ~100-150)
    /// let io_like = ThermalState::calculate(0.015, 4.5, 150.0);
    /// assert!(io_like.tidal_power > io_like.radiogenic_power);
    /// ```
    pub fn calculate(mass_earth: f64, age_gyr: f64, tidal_heating_factor: f64) -> Self {
        let radiogenic_power = Self::calculate_radiogenic_heating(mass_earth, age_gyr);
        let tidal_power = Self::calculate_tidal_heating(mass_earth, tidal_heating_factor);

        let core_temperature =
            Self::estimate_core_temperature(mass_earth, age_gyr, radiogenic_power, tidal_power);

        let mantle_temperature = core_temperature * 0.3;

        let heat_flow =
            Self::estimate_heat_flow(mass_earth, age_gyr, radiogenic_power, tidal_power);

        let heat_sources = Self::determine_heat_sources(age_gyr, radiogenic_power, tidal_power);

        Self {
            core_temperature,
            mantle_temperature,
            heat_flow,
            heat_sources,
            radiogenic_power,
            tidal_power,
        }
    }

    /// Calculate radiogenic heating from long-lived isotopes (U-238, Th-232, K-40)
    fn calculate_radiogenic_heating(mass_earth: f64, age_gyr: f64) -> f64 {
        // Earth's current radiogenic heating: ~8 TW
        const EARTH_RADIOGENIC_PRESENT: f64 = 8.0e12; // W
        const DECAY_TIMESCALE: f64 = 2.5; // Effective half-life in Gyr

        let decay_factor = 2.0_f64.powf(-age_gyr / DECAY_TIMESCALE);
        EARTH_RADIOGENIC_PRESENT * mass_earth * decay_factor
    }

    /// Calculate tidal heating from normalized factor
    fn calculate_tidal_heating(mass_earth: f64, tidal_heating_factor: f64) -> f64 {
        // Baseline calibrated to Io: 1.0e14 W at 0.015 M⊕ with factor ~150
        const BASELINE_POWER_DENSITY: f64 = 4.5e13; // W per Earth mass per unit factor
        BASELINE_POWER_DENSITY * mass_earth * tidal_heating_factor
    }

    /// Estimate core temperature from heat sources
    fn estimate_core_temperature(
        mass_earth: f64,
        age_gyr: f64,
        radiogenic_power: f64,
        tidal_power: f64,
    ) -> f64 {
        // Base temperature from primordial heat
        let base_temp = if mass_earth < 10.0 {
            5000.0 + (mass_earth - 1.0) * 800.0
        } else {
            15000.0 + (mass_earth - 10.0) * 100.0
        };

        // Cooling factor
        let cooling_timescale = mass_earth.powf(0.6) * 3.0;
        let cooling_factor = (-age_gyr / cooling_timescale).exp();

        let min_retention = 0.6;
        let cooled_temp = base_temp * (min_retention + (1.0 - min_retention) * cooling_factor);

        // Additional heating from active sources
        let total_heating = radiogenic_power + tidal_power;
        let heating_contribution = (total_heating / 1.0e12).sqrt() * 500.0;

        cooled_temp + heating_contribution
    }

    /// Estimate surface heat flow
    fn estimate_heat_flow(
        mass_earth: f64,
        age_gyr: f64,
        radiogenic_power: f64,
        tidal_power: f64,
    ) -> f64 {
        // Calculate surface area
        const EARTH_RADIUS: f64 = 6.371e6; // m
        let radius_estimate = EARTH_RADIUS * mass_earth.powf(0.4);
        let surface_area = 4.0 * PI * radius_estimate.powi(2);

        // Primordial heat loss
        let primordial_power = if mass_earth < 10.0 {
            let cooling_timescale = 1.5 + mass_earth.powf(0.8) * 2.0;
            let earth_primordial_ref = 1.67e14;
            let cooling_rate = (-age_gyr / cooling_timescale).exp();
            earth_primordial_ref * mass_earth * cooling_rate
        } else {
            5.0e13 * mass_earth * (-age_gyr / 10.0).exp()
        };

        let total_power = radiogenic_power + tidal_power + primordial_power;
        total_power / surface_area
    }

    /// Determine primary heat sources
    fn determine_heat_sources(
        age_gyr: f64,
        radiogenic_power: f64,
        tidal_power: f64,
    ) -> Vec<HeatSource> {
        let mut sources = Vec::new();

        if age_gyr < 0.1 {
            sources.push(HeatSource::Accretion);
        }

        if tidal_power > radiogenic_power * 2.0 {
            sources.push(HeatSource::Tidal);
            sources.push(HeatSource::Radiogenic);
        } else {
            sources.push(HeatSource::Radiogenic);
            if tidal_power > radiogenic_power * 0.1 {
                sources.push(HeatSource::Tidal);
            }
        }

        sources.push(HeatSource::Primordial);
        sources
    }

    /// Check if tidal heating dominates
    pub fn is_tidally_heated(&self) -> bool {
        self.tidal_power > self.radiogenic_power
    }

    /// Check if planet is geologically active based on heat flow
    pub fn is_geologically_active(&self) -> bool {
        // Earth: ~0.087 W/m², Mars (inactive): ~0.016 W/m²
        self.heat_flow > 0.03
    }
}

impl DifferentiationState {
    /// Determine differentiation state from mass and age
    pub fn classify(mass_earth: f64, age_gyr: f64) -> Self {
        if mass_earth < 0.1 {
            Self::Undifferentiated
        } else if age_gyr < 0.1 {
            Self::PartiallyDifferentiated
        } else {
            Self::FullyDifferentiated
        }
    }
}
