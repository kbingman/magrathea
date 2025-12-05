//! Unified star system output types
//!
//! This crate defines the common output format for all system generation approaches:
//! statistical sampling, stellar-forge formation simulation, and manual construction.

pub mod architecture;
pub mod metadata;
pub mod system;

// Re-export main types at crate root
pub use architecture::SystemArchitecture;
pub use metadata::{GenerationMethod, SystemMetadata};
pub use system::{HabitableZone, PlanetarySystem, snow_line};

// Re-export planetary types for convenience
pub use planetary::composition::Composition;
pub use planetary::planet::{HostStar, Planet};
pub use planetary::planet_class::PlanetClass;
pub use planetary::planet_type::PlanetType;

#[cfg(test)]
mod architecture_test;
#[cfg(test)]
mod metadata_test;
#[cfg(test)]
mod system_test;
