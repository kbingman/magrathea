//! Discrete bodies in protoplanetary disk simulations.
//!
//! This module provides discrete body dynamics for planetesimals, embryos,
//! and protoplanets. Bodies transition from statistical particle populations
//! to individually tracked objects when they become massive enough or when
//! gravitational interactions become important.
//!
//! # Physical Representation
//!
//! Bodies are represented in Cartesian coordinates for numerical integration
//! (compatible with the nbody crate), with orbital elements stored and updated
//! for physics calculations and display.
//!
//! # Integration with nbody
//!
//! DiscreteBody implements the `Massive` trait from the nbody crate, allowing
//! it to work directly with Barnes-Hut tree gravity and symplectic integrators.

mod discrete_body;
mod envelope;
mod orbital_elements;

#[cfg(test)]
mod discrete_body_test;
#[cfg(test)]
mod growth_integration_test;
#[cfg(test)]
mod orbital_elements_test;

pub use discrete_body::{DiscreteBody, EnvelopeState, MassiveBody};
pub use envelope::{
    bondi_radius, can_capture_envelope, critical_core_mass, gas_sound_speed,
    kelvin_helmholtz_growth_rate, kelvin_helmholtz_timescale, supply_limited_accretion_rate,
};
pub use orbital_elements::{
    OrbitalElements, cartesian_to_orbital_elements, orbital_elements_to_cartesian,
};
