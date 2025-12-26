//! Moon (satellite) classification and characterization
//!
//! Moons are planets orbiting planets. They use the same classification system
//! as standalone planets, with additional considerations for tidal heating
//! from their host planet.

use serde::{Deserialize, Serialize};
use units::{Length, Mass};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::planet_class::PlanetClass;
use crate::planet_type::PlanetType;

/// Level of tidal heating from host planet
///
/// Tidal heating arises from gravitational flexing as a moon's orbit
/// is perturbed (often by resonances with other moons).
///
/// # References
/// - Peale et al. (1979) - "Melting of Io by tidal dissipation"
/// - Hussmann et al. (2002) - "Thermal-orbital evolution of Io and Europa"
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum TidalHeatingLevel {
    /// Negligible tidal heating (< 0.1 W/m²)
    Negligible,

    /// Moderate tidal heating (0.1-1.0 W/m², Europa-like)
    /// Subsurface ocean possible, cryovolcanism
    Moderate,

    /// Extreme tidal heating (> 1.0 W/m², Io-like)
    /// Continuous volcanism, surface resurfacing
    Extreme,
}

impl TidalHeatingLevel {
    /// Classify tidal heating level from heat flux in W/m²
    pub fn from_heat_flux(flux: f64) -> Self {
        match flux {
            f if f > 1.0 => Self::Extreme,
            f if f > 0.1 => Self::Moderate,
            _ => Self::Negligible,
        }
    }
}

impl std::fmt::Display for TidalHeatingLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Negligible => write!(f, "Negligible"),
            Self::Moderate => write!(f, "Moderate (Europa-like)"),
            Self::Extreme => write!(f, "Extreme (Io-like)"),
        }
    }
}

/// Formation mechanism for a moon
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum MoonFormation {
    /// Giant impact between planet and large body (Earth-Moon style)
    GiantImpact,

    /// Gravitational capture of passing body (Triton style)
    Capture {
        /// Retrograde orbit (opposite to planet rotation)
        retrograde: bool,
    },

    /// Co-accretion in circumplanetary disk (Galilean style)
    CoAccretion {
        /// Mean-motion resonance with other moons, if any (e.g., 2:1, 4:2:1)
        resonance_order: Option<(u8, u8)>,
    },
}

impl std::fmt::Display for MoonFormation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GiantImpact => write!(f, "Giant Impact"),
            Self::Capture { retrograde: true } => write!(f, "Captured (retrograde)"),
            Self::Capture { retrograde: false } => write!(f, "Captured"),
            Self::CoAccretion {
                resonance_order: Some((a, b)),
            } => write!(f, "Co-accretion ({}:{} resonance)", a, b),
            Self::CoAccretion {
                resonance_order: None,
            } => write!(f, "Co-accretion"),
        }
    }
}

/// A major moon (satellite) with full classification
///
/// Moons are treated as planets orbiting planets. The same `PlanetClass` and
/// `PlanetType` classification applies, with tidal heating factored into
/// temperature and geological activity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct Moon {
    /// Unique identifier within system (e.g., "uuid-b-I")
    pub id: String,

    /// Display name (e.g., "KV-4729 b I")
    pub name: String,

    /// Moon mass
    pub mass: Mass,

    /// Moon radius
    pub radius: Length,

    /// Orbital semi-major axis from host planet center
    pub semi_major_axis: Length,

    /// Orbital eccentricity
    pub eccentricity: f64,

    /// Physical classification (same as planets)
    pub class: PlanetClass,

    /// Observable expression (same as planets)
    pub moon_type: PlanetType,

    /// Formation mechanism
    pub formation: MoonFormation,

    /// Surface/equilibrium temperature in Kelvin (includes tidal + stellar heating)
    pub surface_temp: f64,

    /// Tidal heating flux in W/m²
    pub tidal_heat_flux: f64,

    /// Is tidally locked to host planet?
    pub tidally_locked: bool,
}

impl Moon {
    /// Create a new moon with all properties specified
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: String,
        name: String,
        mass: Mass,
        radius: Length,
        semi_major_axis: Length,
        eccentricity: f64,
        class: PlanetClass,
        moon_type: PlanetType,
        formation: MoonFormation,
        surface_temp: f64,
        tidal_heat_flux: f64,
        tidally_locked: bool,
    ) -> Self {
        Self {
            id,
            name,
            mass,
            radius,
            semi_major_axis,
            eccentricity,
            class,
            moon_type,
            formation,
            surface_temp,
            tidal_heat_flux,
            tidally_locked,
        }
    }

    /// Returns the tidal heating level for this moon
    pub fn tidal_heating_level(&self) -> TidalHeatingLevel {
        TidalHeatingLevel::from_heat_flux(self.tidal_heat_flux)
    }
}

impl std::fmt::Display for Moon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} ({:.3} M⊕, {})",
            self.name,
            self.mass.to_earth_masses(),
            self.moon_type
        )
    }
}

/// A system of moons orbiting a planet
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct MoonSystem {
    /// Major moons (large enough for classification, ~100+ km diameter)
    pub moons: Vec<Moon>,

    /// Has ring system?
    pub has_rings: bool,
}

impl MoonSystem {
    /// Create a new empty moon system
    pub fn new() -> Self {
        Self {
            moons: Vec::new(),
            has_rings: false,
        }
    }

    /// Create a moon system with rings but no major moons
    pub fn rings_only() -> Self {
        Self {
            moons: Vec::new(),
            has_rings: true,
        }
    }

    /// Create a moon system with moons and optional rings
    pub fn with_moons(moons: Vec<Moon>, has_rings: bool) -> Self {
        Self { moons, has_rings }
    }

    /// Add a moon to the system
    pub fn add_moon(&mut self, moon: Moon) {
        self.moons.push(moon);
    }

    /// Returns the number of major moons
    pub fn moon_count(&self) -> usize {
        self.moons.len()
    }

    /// Returns whether this system has any major moons
    pub fn has_moons(&self) -> bool {
        !self.moons.is_empty()
    }

    /// Returns total mass of all moons
    pub fn total_moon_mass(&self) -> Mass {
        self.moons
            .iter()
            .fold(Mass::from_earth_masses(0.0), |acc, m| acc + m.mass)
    }

    /// Sort moons by semi-major axis (innermost first)
    pub fn sort_by_distance(&mut self) {
        self.moons.sort_by(|a, b| {
            a.semi_major_axis
                .to_m()
                .partial_cmp(&b.semi_major_axis.to_m())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}

/// Convert moon index to Roman numeral designation
///
/// Following convention for satellite naming: I, II, III, IV, etc.
/// Ordered by semi-major axis (innermost = I).
///
/// # Examples
/// - 0 → "I"
/// - 1 → "II"
/// - 3 → "IV"
/// - 8 → "IX"
pub fn moon_numeral(index: usize) -> String {
    const NUMERALS: [(usize, &str); 13] = [
        (1000, "M"),
        (900, "CM"),
        (500, "D"),
        (400, "CD"),
        (100, "C"),
        (90, "XC"),
        (50, "L"),
        (40, "XL"),
        (10, "X"),
        (9, "IX"),
        (5, "V"),
        (4, "IV"),
        (1, "I"),
    ];

    NUMERALS
        .iter()
        .fold(
            (String::new(), index + 1),
            |(acc, remaining), &(value, numeral)| {
                let count = remaining / value;
                (acc + &numeral.repeat(count), remaining % value)
            },
        )
        .0
}
