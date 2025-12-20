//! Binary star orbital parameters and stability calculations.
//!
//! This module provides types for representing binary star orbits and calculating
//! stable planetary orbital zones in binary systems.
//!
//! # References
//! - Holman & Wiegert (1999) - "Dynamical Stability of Planets in Binary Systems"
//! - Raghavan et al. (2010) - "A Survey of Stellar Families"

pub mod generation;
pub mod orbital_params;
pub mod stability;

pub use generation::{
    binary_fraction, generate_binary_system, generate_companion, sample_binary_eccentricity,
    sample_mass_ratio, sample_separation,
};
pub use orbital_params::{BinaryConfiguration, BinaryOrbitType, OrbitalParameters};
pub use stability::{
    habitable_zone, p_type_stability_limit, s_type_stability_limit, stable_orbital_range,
};
