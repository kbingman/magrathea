//! Planetary classification and characterization
//!
//! This crate provides types and functions for classifying and characterizing
//! exoplanets based on their physical properties (mass, radius, temperature, etc.).

pub mod composition;
pub mod planet;
pub mod planet_class;
pub mod planet_type;

// Re-export key types at crate root
pub use composition::Composition;
pub use planet::{HostStar, Planet};
pub use planet_class::PlanetClass;
pub use planet_type::PlanetType;

#[cfg(test)]
mod composition_test;
#[cfg(test)]
mod planet_class_test;
#[cfg(test)]
mod planet_test;
#[cfg(test)]
mod planet_type_test;
