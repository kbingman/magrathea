//! Secondary planet classification: observable expressions
//!
//! PlanetType represents how a planet of a given PlanetClass manifests
//! based on temperature, composition, stellar environment, and history.
//! This is the second tier of the two-tier classification system.

use serde::{Deserialize, Serialize};

use crate::composition::Composition;
use crate::planet_class::PlanetClass;

/// Temperature thresholds for type determination (Kelvin)
pub mod temperature {
    /// Above this: lava world with magma ocean
    pub const LAVA: f64 = 1500.0;
    /// Above this: hot Jupiter regime
    pub const HOT_JUPITER: f64 = 1000.0;
    /// Above this: warm, no surface ice
    pub const WARM: f64 = 350.0;
    /// Below this: frozen world
    pub const FROZEN: f64 = 200.0;
    /// Habitable zone bounds (approximate)
    pub const HZ_INNER: f64 = 350.0;
    pub const HZ_OUTER: f64 = 200.0;
}

/// Flux thresholds (in Earth flux units, F⊕)
pub mod flux {
    /// Above this, transitional planets likely stripped to rocky core
    pub const EVAPORATION_THRESHOLD: f64 = 100.0;
}

/// Estimate whether a planet is tidally locked based on orbital period and stellar mass
///
/// Tidal locking timescale scales as τ ∝ a^6/M_star. For planets in the habitable zone,
/// the threshold period varies by stellar type because the HZ location scales with
/// luminosity while locking timescale scales differently.
///
/// Thresholds are calibrated so that:
/// - M-dwarf HZ planets (P ~ 10-40 days) are typically locked
/// - K-dwarf HZ planets (P ~ 50-100 days) are sometimes locked
/// - G/F-dwarf HZ planets (P ~ 200-400 days) are rarely locked
fn is_tidally_locked(orbital_period_days: f64, stellar_mass_solar: f64) -> bool {
    let threshold_days = match stellar_mass_solar {
        m if m < 0.20 => 20.0, // Late M: very close HZ
        m if m < 0.35 => 30.0, // Mid M
        m if m < 0.50 => 40.0, // Early M
        m if m < 0.70 => 15.0, // K-dwarf: HZ further out, harder to lock
        _ => 5.0,              // G/F: almost never locked in HZ
    };
    orbital_period_days < threshold_days
}

/// Observable planet type - the expression of a physical regime
///
/// Each variant represents a distinct observable category determined by
/// the planet's mass regime (PlanetClass) combined with its temperature,
/// composition, and stellar environment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "camelCase")]
pub enum PlanetType {
    // =========================================================================
    // Rocky expressions (PlanetClass::Rocky, M < 2 M⊕)
    // =========================================================================
    /// M < 0.5 M⊕: Mars/Mercury class, cannot retain atmosphere long-term
    SubEarth,

    /// Airless, cratered surface (Mercury, Moon-like)
    Barren,

    /// T > 1500 K: Permanent magma ocean, silicate/metal vapor atmosphere
    Lava {
        tidally_locked: bool,
        surface_temp_k: f64,
    },

    /// 300-500 K: Thin CO₂ or N₂ atmosphere (Mars/Venus-like)
    Desert,

    /// T < 200 K: Ice-covered, possible subsurface ocean
    Frozen { has_subsurface_ocean: bool },

    /// High water fraction, surface liquid water ocean
    Oceanic { ocean_fraction: f64 },

    /// Earth-like: temperate, partial ocean coverage, active geology
    Terran {
        ocean_fraction: f64,
        tectonically_active: bool,
    },

    /// Tidally locked with habitable terminator zone
    Eyeball { terminator_habitable: bool },

    /// C/O > 0.8: Carbon-dominated surface (graphite, carbides, diamond)
    Carbon,

    /// Exposed metallic core (>70% iron, mantle-stripped)
    Iron,

    // =========================================================================
    // Transitional expressions (PlanetClass::Transitional, 2-5 M⊕)
    // =========================================================================
    /// High stellar flux or stripped: massive rocky world, secondary atmosphere
    SuperTerran,

    /// Low stellar flux, retained H/He envelope (1-10% by mass)
    MiniNeptune { envelope_fraction: f64 },

    /// High volatile fraction: deep global ocean, possibly over high-P ice
    WaterWorld { has_hp_ice_layer: bool },

    /// H₂ atmosphere over liquid water ocean - astrobiologically interesting
    Hycean,

    /// Stripped gas giant core (M > 10 M⊕ rocky, remnant of evaporated giant)
    Chthonian,

    // =========================================================================
    // Volatile expressions (PlanetClass::Volatile, 5-160 M⊕)
    // =========================================================================
    /// Warm, close-in, puffy hydrogen-dominated
    SubNeptune,

    /// Intermediate distance, transitional cloud layers
    WarmNeptune,

    /// Cold outer system: Uranus/Neptune-like, methane clouds
    IceGiant { has_rings: bool },

    // =========================================================================
    // Giant expressions (PlanetClass::Giant, > 160 M⊕)
    // =========================================================================
    /// Cold Jupiter-like: banded clouds, internal heat source
    GasGiant { has_rings: bool },

    /// T > 1000 K: Inflated, tidally locked, atmospheric escape
    HotJupiter {
        equilibrium_temp_k: f64,
        evaporating: bool,
    },

    /// Inflated sub-Jovian with anomalously low density (0.1-0.5 M_J)
    PuffySaturn,
}

impl PlanetType {
    /// Returns the parent PlanetClass for this type
    pub fn class(&self) -> PlanetClass {
        match self {
            // Rocky expressions
            Self::SubEarth
            | Self::Barren
            | Self::Lava { .. }
            | Self::Desert
            | Self::Frozen { .. }
            | Self::Oceanic { .. }
            | Self::Terran { .. }
            | Self::Eyeball { .. }
            | Self::Carbon
            | Self::Iron => PlanetClass::Rocky,

            // Transitional expressions
            Self::SuperTerran
            | Self::MiniNeptune { .. }
            | Self::WaterWorld { .. }
            | Self::Hycean
            | Self::Chthonian => PlanetClass::Transitional,

            // Volatile expressions
            Self::SubNeptune | Self::WarmNeptune | Self::IceGiant { .. } => PlanetClass::Volatile,

            // Giant expressions
            Self::GasGiant { .. } | Self::HotJupiter { .. } | Self::PuffySaturn => {
                PlanetClass::Giant
            }
        }
    }

    /// Determine planet type from class, composition, and environment
    ///
    /// This is the core mapping function that takes a physical classification
    /// and observable conditions to determine the planet's expression.
    ///
    /// # Arguments
    /// * `class` - Physical mass regime
    /// * `composition` - Bulk composition fractions
    /// * `equilibrium_temp` - Equilibrium temperature in Kelvin
    /// * `incident_flux` - Stellar flux in Earth flux units (F⊕)
    /// * `mass_earth` - Planet mass in Earth masses
    /// * `stellar_mass` - Stellar mass in solar masses (for tidal locking estimate)
    /// * `semi_major_axis_au` - Orbital distance in AU (for tidal locking estimate)
    pub fn from_environment(
        class: PlanetClass,
        composition: &Composition,
        equilibrium_temp: f64,
        incident_flux: f64,
        mass_earth: f64,
        stellar_mass: f64,
        semi_major_axis_au: f64,
    ) -> Self {
        match class {
            PlanetClass::Rocky => Self::determine_rocky(
                composition,
                equilibrium_temp,
                mass_earth,
                incident_flux,
                stellar_mass,
                semi_major_axis_au,
            ),
            PlanetClass::Transitional => {
                Self::determine_transitional(composition, equilibrium_temp, incident_flux)
            }
            PlanetClass::Volatile => Self::determine_volatile(equilibrium_temp, mass_earth),
            PlanetClass::Giant => Self::determine_giant(equilibrium_temp, mass_earth),
        }
    }

    fn determine_rocky(
        composition: &Composition,
        t_eq: f64,
        mass_earth: f64,
        _incident_flux: f64,
        stellar_mass: f64,
        semi_major_axis_au: f64,
    ) -> Self {
        // Check for frozen worlds first - TNOs and other icy bodies
        // This takes precedence over mass check since frozen dwarf planets
        // (Pluto, Eris, etc.) are distinct from warm sub-Earths
        if t_eq < temperature::FROZEN {
            return Self::Frozen {
                has_subsurface_ocean: composition.water > 0.05,
            };
        }

        // Check mass - very small warm bodies
        if mass_earth < 0.5 {
            return Self::SubEarth;
        }

        // Check extreme compositions
        if composition.iron > 0.6 {
            return Self::Iron;
        }

        // Temperature-based determination
        match t_eq {
            t if t > temperature::LAVA => Self::Lava {
                tidally_locked: true,
                surface_temp_k: t,
            },
            t if t < temperature::FROZEN => Self::Frozen {
                has_subsurface_ocean: composition.water > 0.05,
            },
            t if t > temperature::HZ_INNER => Self::Desert,
            t if t >= temperature::HZ_OUTER => {
                // Habitable zone - check for tidal locking (Eyeball worlds)
                // Calculate orbital period: P = sqrt(a³/M_star) in years, convert to days
                let orbital_period_days =
                    (semi_major_axis_au.powi(3) / stellar_mass).sqrt() * 365.25;
                let likely_tidally_locked = is_tidally_locked(orbital_period_days, stellar_mass);

                if likely_tidally_locked && composition.water > 0.001 {
                    // Eyeball world: tidally locked with potential terminator habitability
                    Self::Eyeball {
                        terminator_habitable: composition.water > 0.01,
                    }
                } else if composition.water > 0.3 {
                    Self::Oceanic {
                        ocean_fraction: 0.95,
                    }
                } else if composition.water > 0.001 {
                    Self::Terran {
                        ocean_fraction: 0.7,
                        tectonically_active: mass_earth > 0.5,
                    }
                } else {
                    Self::Desert
                }
            }
            _ => Self::Frozen {
                has_subsurface_ocean: false,
            },
        }
    }

    fn determine_transitional(composition: &Composition, t_eq: f64, incident_flux: f64) -> Self {
        // Photoevaporation determines envelope fate
        let stripped = incident_flux > flux::EVAPORATION_THRESHOLD;

        if stripped || composition.h_he_gas <= 0.0 {
            // Stripped or never had envelope
            if composition.water > 0.3 {
                Self::WaterWorld {
                    has_hp_ice_layer: composition.water > 0.5,
                }
            } else {
                Self::SuperTerran
            }
        } else {
            // Retained envelope
            if composition.water > 0.3
                && t_eq > temperature::HZ_OUTER
                && t_eq < temperature::HZ_INNER
            {
                Self::Hycean
            } else {
                Self::MiniNeptune {
                    envelope_fraction: composition.h_he_gas,
                }
            }
        }
    }

    fn determine_volatile(t_eq: f64, mass_earth: f64) -> Self {
        // Saturn-class (50-160 M⊕) with moderate irradiation → Puffy Saturn
        // These are inflated sub-Jovians with anomalously low density
        let is_saturn_mass = mass_earth > 50.0;
        let is_warm_enough_for_inflation = t_eq > 300.0 && t_eq < 1000.0;

        if is_saturn_mass && is_warm_enough_for_inflation {
            return Self::PuffySaturn;
        }

        match t_eq {
            t if t > 400.0 => Self::SubNeptune,
            t if t > 150.0 => Self::WarmNeptune,
            _ => Self::IceGiant { has_rings: false },
        }
    }

    fn determine_giant(t_eq: f64, mass_earth: f64) -> Self {
        let mass_jupiter = mass_earth / 317.8;

        match (t_eq, mass_jupiter) {
            (t, _) if t > temperature::HOT_JUPITER => Self::HotJupiter {
                equilibrium_temp_k: t,
                evaporating: t > 1500.0,
            },
            (_, m) if m < 0.5 => Self::PuffySaturn,
            _ => Self::GasGiant { has_rings: false },
        }
    }

    /// Returns whether this type typically has a solid surface
    pub fn has_solid_surface(&self) -> bool {
        match self {
            Self::SubEarth
            | Self::Barren
            | Self::Lava { .. }
            | Self::Desert
            | Self::Frozen { .. }
            | Self::Terran { .. }
            | Self::Eyeball { .. }
            | Self::Carbon
            | Self::Iron
            | Self::SuperTerran
            | Self::Chthonian => true,

            Self::Oceanic { .. } | Self::WaterWorld { .. } | Self::Hycean => false, // Surface is liquid

            Self::MiniNeptune { .. }
            | Self::SubNeptune
            | Self::WarmNeptune
            | Self::IceGiant { .. }
            | Self::GasGiant { .. }
            | Self::HotJupiter { .. }
            | Self::PuffySaturn => false,
        }
    }

    /// Returns whether this type could potentially support surface liquid water
    pub fn potentially_habitable(&self) -> bool {
        matches!(
            self,
            Self::Terran { .. }
                | Self::Oceanic { .. }
                | Self::Eyeball {
                    terminator_habitable: true
                }
                | Self::Hycean
        )
    }

    /// Human-readable description
    pub fn description(&self) -> &'static str {
        match self {
            Self::SubEarth => "Sub-Earth: small rocky body, cannot retain atmosphere",
            Self::Barren => "Barren: airless, cratered rocky world",
            Self::Lava { .. } => "Lava world: permanent magma ocean",
            Self::Desert => "Desert: thin atmosphere, no surface liquid",
            Self::Frozen { .. } => "Frozen: ice-covered, possible subsurface ocean",
            Self::Oceanic { .. } => "Oceanic: water-dominated surface",
            Self::Terran { .. } => "Terran: Earth-like with oceans and continents",
            Self::Eyeball { .. } => "Eyeball: tidally locked, terminator habitability",
            Self::Carbon => "Carbon world: graphite/carbide surface",
            Self::Iron => "Iron world: exposed metallic core",
            Self::SuperTerran => "Super-Terran: massive rocky world",
            Self::MiniNeptune { .. } => "Mini-Neptune: rocky core with H/He envelope",
            Self::WaterWorld { .. } => "Water world: deep global ocean",
            Self::Hycean => "Hycean: ocean under hydrogen atmosphere",
            Self::Chthonian => "Chthonian: stripped gas giant core",
            Self::SubNeptune => "Sub-Neptune: warm, puffy, hydrogen-rich",
            Self::WarmNeptune => "Warm Neptune: intermediate temperature",
            Self::IceGiant { .. } => "Ice giant: cold, methane clouds",
            Self::GasGiant { .. } => "Gas giant: Jupiter-like, banded clouds",
            Self::HotJupiter { .. } => "Hot Jupiter: close-in, inflated",
            Self::PuffySaturn => "Puffy Saturn: low-density sub-Jovian",
        }
    }
}

impl std::fmt::Display for PlanetType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::SubEarth => "Sub-Earth",
            Self::Barren => "Barren",
            Self::Lava { .. } => "Lava World",
            Self::Desert => "Desert",
            Self::Frozen { .. } => "Frozen",
            Self::Oceanic { .. } => "Oceanic",
            Self::Terran { .. } => "Terran",
            Self::Eyeball { .. } => "Eyeball",
            Self::Carbon => "Carbon World",
            Self::Iron => "Iron World",
            Self::SuperTerran => "Super-Terran",
            Self::MiniNeptune { .. } => "Mini-Neptune",
            Self::WaterWorld { .. } => "Water World",
            Self::Hycean => "Hycean",
            Self::Chthonian => "Chthonian",
            Self::SubNeptune => "Sub-Neptune",
            Self::WarmNeptune => "Warm Neptune",
            Self::IceGiant { .. } => "Ice Giant",
            Self::GasGiant { .. } => "Gas Giant",
            Self::HotJupiter { .. } => "Hot Jupiter",
            Self::PuffySaturn => "Puffy Saturn",
        };
        write!(f, "{}", name)
    }
}
