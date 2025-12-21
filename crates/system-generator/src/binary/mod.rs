//! Binary star generation for stellar systems.
//!
//! This module provides functions for generating binary star companions
//! and determining binary fractions based on stellar mass.
//!
//! # References
//! - Raghavan et al. (2010) - "A Survey of Stellar Families"

pub mod generation;

pub use generation::{
    binary_fraction, generate_binary_system, generate_companion, sample_binary_eccentricity,
    sample_mass_ratio, sample_separation,
};

// Re-export binary types from star-system for convenience
pub use star_system::binary::{
    BinaryConfiguration, BinaryOrbitType, OrbitalParameters, habitable_zone,
    p_type_stability_limit, s_type_stability_limit, stable_orbital_range,
};
