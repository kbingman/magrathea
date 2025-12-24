//! Planetary classification and characterization
//!
//! This crate provides types and functions for classifying and characterizing
//! exoplanets based on their physical properties (mass, radius, temperature, etc.).
//!
//! # Two-Tier Classification System
//!
//! Planets are classified using two complementary systems:
//!
//! 1. **[`PlanetClass`]** - Physical mass regime (Compact, Transitional, Volatile, Giant)
//! 2. **[`PlanetType`]** - Observable expression (21+ variants based on environment)
//!
//! # Additional Characterization
//!
//! - **[`AtmosphereType`]** - Atmospheric composition classification
//! - **[`Temperature`]** - Equilibrium and effective temperature calculations
//! - **[`Composition`]** - Bulk composition (iron, rock, water, gas fractions)

pub mod atmosphere;
pub mod composition;
pub mod planet;
pub mod planet_class;
pub mod planet_type;
pub mod temperature;

// Re-export key types at crate root
pub use atmosphere::AtmosphereType;
pub use composition::Composition;
pub use planet::{HostStar, Planet};
pub use planet_class::PlanetClass;
pub use planet_type::PlanetType;
pub use temperature::{Temperature, TemperatureClass};

#[cfg(test)]
mod atmosphere_test;
#[cfg(test)]
mod composition_test;
#[cfg(test)]
mod planet_class_test;
#[cfg(test)]
mod planet_test;
#[cfg(test)]
mod planet_type_test;
#[cfg(test)]
mod temperature_test;
