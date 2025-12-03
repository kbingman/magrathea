//! Stellar object types and generation functions.

pub mod generation;
pub mod sampling;
pub mod spectral;
pub mod stellar_color;
pub mod stellar_objects;
pub mod stellar_radius;

#[cfg(test)]
mod generation_test;
#[cfg(test)]
mod sampling_test;
#[cfg(test)]
mod spectral_test;
#[cfg(test)]
mod stellar_color_test;

// Re-export types
pub use stellar_objects::{
    BlackHole, GiantStar, MainSequenceStar, NeutronStar, StellarObject, WhiteDwarf, WhiteDwarfType,
};
pub use stellar_radius::StellarRadius;

// Re-export spectral types
pub use spectral::{LuminosityClass, SpectralType, VariabilityType};

// Re-export color
pub use stellar_color::StellarColor;

// Re-export generation functions
pub use generation::{
    black_hole, giant_star, main_sequence_star, neutron_star, sample_main_sequence_star,
    sample_stellar_object, solar_analog, stellar_object, white_dwarf,
};
