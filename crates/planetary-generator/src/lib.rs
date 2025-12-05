//! Statistical planetary system generation
//!
//! Generates statistically realistic planetary systems based on occurrence rates
//! from Kepler/TESS, radial velocity, microlensing, and direct imaging surveys.

pub mod generation;
pub mod sampling;

// Re-export main generation functions
pub use generation::{
    from_star, from_star_with_id, generate_planetary_system, generate_planetary_system_named,
    generate_planetary_system_random,
};

// Re-export sampling functions
pub use sampling::{
    period_to_semi_major_axis, sample_eccentricity, sample_inclination, sample_orbital_period,
    sample_planet_count, sample_planet_mass,
};

// Re-export star-system types for convenience
pub use star_system::{
    GenerationMethod, HabitableZone, PlanetarySystem, SystemArchitecture, SystemMetadata,
};

#[cfg(test)]
mod generation_test;
#[cfg(test)]
mod sampling_test;
