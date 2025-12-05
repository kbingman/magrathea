use std::fmt;

use serde::{Deserialize, Serialize};

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum SpectralType {
    // Main sequence types
    O,
    B,
    A,
    F,
    G,
    K,
    M,
    L,
    T,
    Y,  // Brown dwarfs
    D,  // White dwarf
    N,  // Neutron star
    BH, // Black hole

    // Other rare types
    S, // S-type (similar to M but with zirconium oxide bands)
    C, // Carbon stars
    R, // Late-type carbon stars
    Q, // Nova remnants
}

impl fmt::Display for SpectralType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let str = match self {
            SpectralType::O => "O",
            SpectralType::B => "B",
            SpectralType::A => "A",
            SpectralType::F => "F",
            SpectralType::G => "G",
            SpectralType::K => "K",
            SpectralType::M => "M",
            SpectralType::L => "L",
            SpectralType::T => "T",
            SpectralType::Y => "Y",
            SpectralType::D => "D",
            SpectralType::N => "N",
            SpectralType::BH => "BH",
            SpectralType::S => "S",
            SpectralType::C => "C",
            SpectralType::R => "R",
            SpectralType::Q => "Q",
        };
        write!(f, "{}", str)
    }
}

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum LuminosityClass {
    IAPLUS, // Hypergiants
    IA,     // Bright supergiants
    IB,     // Supergiants
    II,     // Bright giants
    III,    // Normal giants
    IV,     // Subgiants
    V,      // Main sequence
    VI,     // Subdwarfs
    D,      // White dwarf
    BD,     // Brown dwarf
    NS,     // Neutron star
}

impl fmt::Display for LuminosityClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let str = match self {
            LuminosityClass::IAPLUS => "Ia+",
            LuminosityClass::IA => "Ia",
            LuminosityClass::IB => "Ib",
            LuminosityClass::II => "II",
            LuminosityClass::III => "III",
            LuminosityClass::IV => "IV",
            LuminosityClass::V => "V",
            LuminosityClass::VI => "VI",
            LuminosityClass::D => "D",
            LuminosityClass::BD => "BD",
            LuminosityClass::NS => "NS",
        };
        write!(f, "{}", str)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum VariabilityType {
    None,
    Cepheid,     // Classical pulsating supergiants
    RRLyrae,     // Pulsating horizontal branch stars
    Mira,        // Long-period variables
    SemiRegular, // SR variables
    Irregular,   // Irregular variables
    Beta,        // Beta Cephei variables
    Delta,       // Delta Scuti variables
    Flare,       // Flare stars
    Pulsar,      // Pulsar neutron stars
}

impl VariabilityType {
    pub fn determine_variability(mass: f64, temperature: f64, lum_class: LuminosityClass) -> Self {
        match (mass, temperature, lum_class) {
            // Cepheid variables: supergiants in instability strip
            (m, t, LuminosityClass::IA) | (m, t, LuminosityClass::IB)
                if m > 5.0 && t > 5000.0 && t < 7000.0 =>
            {
                VariabilityType::Cepheid
            }

            // RR Lyrae: horizontal branch stars
            (_, t, _) if t > 6000.0 && t < 7500.0 => VariabilityType::RRLyrae,

            // Mira variables: cool giants
            (_, t, LuminosityClass::III) if t < 3500.0 => VariabilityType::Mira,

            // Beta Cephei: hot massive stars
            (m, t, _) if m > 8.0 && t > 20000.0 => VariabilityType::Beta,

            // Delta Scuti: A-F dwarfs
            (_, t, LuminosityClass::V) if t > 6500.0 && t < 8500.0 => VariabilityType::Delta,

            // Flare stars: typically M dwarfs
            (_, t, LuminosityClass::V) if t < 3500.0 => VariabilityType::Flare,

            // Default to non-variable
            _ => VariabilityType::None,
        }
    }
}
