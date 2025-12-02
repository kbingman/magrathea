//! Particle aerodynamics for protoplanetary disks.
//!
//! This module provides the physics of solid particles in gas disks:
//! - Drag regimes (Epstein, Stokes)
//! - Stopping times and Stokes numbers
//! - Radial drift and vertical settling (future)
//!
//! # Physics
//!
//! Particles in a gas disk experience aerodynamic drag. The drag regime
//! depends on the particle size relative to the gas mean free path:
//!
//! - **Epstein drag**: particle size < 9λ/4 (most disk particles)
//! - **Stokes drag**: particle size > 9λ/4 and Re < 1
//!
//! The stopping time t_s measures how quickly a particle couples to the gas.
//! The Stokes number τ_s = t_s × Ω_K determines aerodynamic behavior:
//! - τ_s << 1: tightly coupled to gas
//! - τ_s ~ 1: maximum radial drift
//! - τ_s >> 1: decoupled from gas

mod particle;

#[cfg(test)]
mod particle_test;

pub use particle::{DragRegime, Particle};
