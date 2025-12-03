//! Planetary composition types
//!
//! Composition is tracked as mass fractions of major components.
//! Field names match existing codebase conventions.

use serde::{Deserialize, Serialize};

/// Bulk composition as mass fractions
///
/// Fractions should sum to approximately 1.0.
/// Uses `h_he_gas` naming to match existing codebase.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Composition {
    /// Iron/nickel core fraction (Earth ~0.32)
    pub iron: f64,
    /// Silicate mantle/crust fraction (Earth ~0.68)
    pub silicate: f64,
    /// Water/ice fraction (Earth ~0.0002, Europa ~0.06)
    pub water: f64,
    /// H/He envelope fraction (Earth ~0, Neptune ~0.15, Jupiter ~0.95)
    /// Named to match existing codebase
    pub h_he_gas: f64,
}

impl Composition {
    /// Create a new composition, normalizing fractions to sum to 1.0
    pub fn new(iron: f64, silicate: f64, water: f64, h_he_gas: f64) -> Self {
        let total = iron + silicate + water + h_he_gas;
        if total <= 0.0 {
            return Self::earth_like();
        }
        Self {
            iron: iron / total,
            silicate: silicate / total,
            water: water / total,
            h_he_gas: h_he_gas / total,
        }
    }

    /// Earth-like composition (rocky, minimal volatiles)
    pub fn earth_like() -> Self {
        Self {
            iron: 0.32,
            silicate: 0.68,
            water: 0.0002,
            h_he_gas: 0.0,
        }
    }

    /// Mercury-like composition (iron-rich)
    pub fn iron_rich() -> Self {
        Self {
            iron: 0.70,
            silicate: 0.30,
            water: 0.0,
            h_he_gas: 0.0,
        }
    }

    /// Mars-like composition
    pub fn mars_like() -> Self {
        Self {
            iron: 0.25,
            silicate: 0.73,
            water: 0.02,
            h_he_gas: 0.0,
        }
    }

    /// Water world composition (high ice fraction)
    pub fn water_world() -> Self {
        Self {
            iron: 0.10,
            silicate: 0.30,
            water: 0.60,
            h_he_gas: 0.0,
        }
    }

    /// Mini-Neptune composition (rocky core with H/He envelope)
    pub fn mini_neptune() -> Self {
        Self {
            iron: 0.15,
            silicate: 0.35,
            water: 0.10,
            h_he_gas: 0.40,
        }
    }

    /// Ice giant composition (Uranus/Neptune-like)
    pub fn ice_giant() -> Self {
        Self {
            iron: 0.02,
            silicate: 0.13,
            water: 0.70,
            h_he_gas: 0.15,
        }
    }

    /// Gas giant composition (Jupiter/Saturn-like)
    pub fn gas_giant() -> Self {
        Self {
            iron: 0.01,
            silicate: 0.02,
            water: 0.02,
            h_he_gas: 0.95,
        }
    }

    /// Sample ice giant composition with stochastic variation
    ///
    /// Ice giants show diversity: Uranus is ~80% ices, Neptune ~65%.
    /// Water/ice fraction varies from 0.50-0.80, H/He from 0.10-0.25.
    pub fn sample_ice_giant(rng: &mut impl rand::Rng) -> Self {
        // Core rock fraction (small)
        let rock = 0.10 + rng.random::<f64>() * 0.10; // 0.10-0.20

        // H/He envelope (modest)
        let envelope = 0.10 + rng.random::<f64>() * 0.15; // 0.10-0.25

        // Water/ice dominates the rest
        let water = 1.0 - rock - envelope;

        // Split rock into iron and silicate
        let iron_frac = 0.15 + rng.random::<f64>() * 0.10; // 15-25% of rock is iron
        Self::new(rock * iron_frac, rock * (1.0 - iron_frac), water, envelope)
    }

    /// Sample gas giant composition with stochastic variation
    ///
    /// Gas giants vary in core mass and metallicity enrichment.
    /// Jupiter: ~3-15 M⊕ core, Saturn: ~15-25 M⊕ core.
    /// H/He fraction ranges from 0.85-0.97.
    pub fn sample_gas_giant(rng: &mut impl rand::Rng) -> Self {
        // H/He dominates
        let envelope = 0.85 + rng.random::<f64>() * 0.12; // 0.85-0.97

        // Remaining mass is core (metals + ices)
        let core = 1.0 - envelope;

        // Core composition varies - some ice, some rock
        let ice_frac = 0.3 + rng.random::<f64>() * 0.4; // 30-70% of core is ice
        let water = core * ice_frac;
        let rock = core * (1.0 - ice_frac);

        // Split rock into iron and silicate
        let iron_frac = 0.20 + rng.random::<f64>() * 0.15;
        Self::new(rock * iron_frac, rock * (1.0 - iron_frac), water, envelope)
    }

    /// Return composition with envelope stripped (photoevaporation)
    pub fn stripped_envelope(&self) -> Self {
        let core_total = self.iron + self.silicate + self.water;
        if core_total <= 0.0 {
            return Self::earth_like();
        }

        Self {
            iron: self.iron / core_total,
            silicate: self.silicate / core_total,
            water: self.water / core_total,
            h_he_gas: 0.0,
        }
    }

    /// Returns the rocky (iron + silicate) fraction
    pub fn rocky_fraction(&self) -> f64 {
        self.iron + self.silicate
    }

    /// Returns the volatile (water + H/He) fraction
    pub fn volatile_fraction(&self) -> f64 {
        self.water + self.h_he_gas
    }

    /// Returns the condensate (everything except H/He) fraction - the "core"
    pub fn core_fraction(&self) -> f64 {
        1.0 - self.h_he_gas
    }

    /// Check if composition is dominated by rock (iron + silicate > 0.5)
    pub fn is_rocky(&self) -> bool {
        self.rocky_fraction() > 0.5
    }

    /// Check if composition is dominated by water (water > 0.25)
    pub fn is_water_rich(&self) -> bool {
        self.water > 0.25
    }

    /// Check if composition has significant H/He envelope (h_he_gas > 0.01)
    pub fn has_envelope(&self) -> bool {
        self.h_he_gas > 0.01
    }

    /// Check if composition is gas-dominated (h_he_gas > 0.5)
    pub fn is_gas_dominated(&self) -> bool {
        self.h_he_gas > 0.5
    }
}

impl Default for Composition {
    fn default() -> Self {
        Self::earth_like()
    }
}
