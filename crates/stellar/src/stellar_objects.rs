//! Pure data types for stellar objects.
//!
//! These structs represent the physical properties of stars and stellar remnants.
//! They are designed as DTOs (data transfer objects) for serialization to JSON
//! and consumption by UI/WASM clients.

use serde::{Deserialize, Serialize};
use units::time::Time;
use units::{Mass, Temperature};

use super::spectral::{LuminosityClass, SpectralType, VariabilityType};
use super::stellar_color::StellarColor;
use super::stellar_radius::StellarRadius;

/// Main sequence stars: core hydrogen burning stars
///
/// These are the most common type of star and the primary hosts for planetary systems.
/// Properties are derived from mass using standard stellar mass-luminosity and
/// mass-temperature relations.
///
/// # Physical Ranges
/// * mass: 0.08-150 solar masses
/// * temperature: 2400K (M-type) to 50,000K (O-type)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MainSequenceStar {
    pub mass: Mass,
    pub radius: StellarRadius,
    /// Luminosity in solar luminosities (L☉)
    pub luminosity: f64,
    pub temperature: Temperature,
    pub spectral_type: SpectralType,
    pub luminosity_class: LuminosityClass,
    /// Spectral subtype (0-9)
    pub subtype: u8,
    pub variability: VariabilityType,
    /// Metallicity [Fe/H] in dex (0.0 = solar)
    pub metallicity: f64,
    pub age: Time,
    /// Visual color based on temperature
    pub color: StellarColor,
}

/// Giant stars: evolved stars burning heavier elements
///
/// Giants have exhausted core hydrogen and expanded significantly.
/// Includes subgiants, giants, bright giants, supergiants, and hypergiants.
///
/// # Physical Ranges
/// * mass: 0.3-8 solar masses (current), from 0.8-100+ solar masses (initial)
/// * radius: 10-1000+ times main sequence radius
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GiantStar {
    pub mass: Mass,
    pub radius: StellarRadius,
    /// Luminosity in solar luminosities (L☉)
    pub luminosity: f64,
    pub temperature: Temperature,
    pub spectral_type: SpectralType,
    /// Spectral subtype (0-9)
    pub subtype: u8,
    pub luminosity_class: LuminosityClass,
    pub variability: VariabilityType,
    /// Visual color based on temperature
    pub color: StellarColor,
}

/// White dwarfs: dense stellar remnants of low/medium mass stars
///
/// The final evolutionary state for stars with initial mass < 8 M☉.
/// Supported by electron degeneracy pressure.
///
/// # Physical Ranges
/// * mass: 0.17-1.44 solar masses (Chandrasekhar limit)
/// * radius: ~0.01 solar radii (Earth-sized)
/// * temperature: 4,000K to 150,000K (cooling over time)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WhiteDwarf {
    pub mass: Mass,
    pub radius: StellarRadius,
    /// Luminosity in solar luminosities (L☉)
    pub luminosity: f64,
    pub temperature: Temperature,
    pub spectral_type: WhiteDwarfType,
    /// Visual color based on temperature
    pub color: StellarColor,
}

/// White dwarf spectral classification
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
pub enum WhiteDwarfType {
    /// Hydrogen-rich atmosphere (most common)
    DA,
    /// Helium-rich atmosphere
    DB,
    /// No strong spectral lines (continuous spectrum)
    DC,
    /// Carbon features
    DQ,
    /// Metal-rich atmosphere
    DZ,
}

/// Neutron stars: ultra-dense stellar remnants
///
/// Formed from core-collapse supernovae of massive stars (8-20 M☉).
/// Supported by neutron degeneracy pressure.
///
/// # Physical Ranges
/// * mass: 1.1-2.5 solar masses
/// * radius: ~10-14 km
/// * magnetic_field: 10^8 to 10^15 Gauss
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NeutronStar {
    pub mass: Mass,
    pub radius: StellarRadius,
    /// Magnetic field strength as log₁₀(B/Gauss)
    pub magnetic_field: f64,
    /// Whether this neutron star is observed as a pulsar
    pub pulsar: bool,
    /// Whether this is a magnetar (B > 10^14 G)
    pub magnetar: bool,
    /// Visual color (blue-white, varies with magnetar/pulsar status)
    pub color: StellarColor,
}

/// Black holes: objects with gravity so strong that nothing can escape
///
/// Formed from the most massive stars (> 20 M☉) or through mergers.
///
/// # Physical Ranges
/// * mass: 3-100 solar masses (stellar), millions to billions (supermassive)
/// * spin: 0 (Schwarzschild) to 1 (maximally rotating Kerr)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BlackHole {
    pub mass: Mass,
    /// Dimensionless spin parameter (0 to 1)
    pub spin: f64,
    /// Whether an accretion disk is present
    pub has_accretion: bool,
    /// Visual color (dark, or orange-white with accretion)
    pub color: StellarColor,
}

/// A stellar object of any type
///
/// This enum wraps all possible stellar objects for use in heterogeneous
/// collections (e.g., a stellar system with mixed object types).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum StellarObject {
    MainSequence(MainSequenceStar),
    Giant(GiantStar),
    WhiteDwarf(WhiteDwarf),
    NeutronStar(NeutronStar),
    BlackHole(BlackHole),
}
