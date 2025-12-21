//! Particle aerodynamics for protoplanetary disks.
//!
//! This module provides the physics of solid particles in gas disks:
//! - Drag regimes (Epstein, Stokes)
//! - Stopping times and Stokes numbers
//! - Radial drift and vertical settling
//! - Size distributions for particle populations
//! - Coagulation (sticking) for dust growth
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
//!
//! # Size Distributions
//!
//! Particle populations are represented statistically via size distributions.
//! The MRN distribution (q = 3.5) is typical for ISM dust. Coagulation
//! evolves the distribution toward larger sizes.
//!
//! # Coagulation
//!
//! Dust growth proceeds through collisions. At low velocities, particles
//! stick to form larger aggregates. The collision rate depends on
//! geometric cross-sections and relative velocities (differential drift,
//! turbulence). The Smoluchowski equation governs the evolution.

mod coagulation;
mod collision_outcome;
mod particle;
mod particle_bin;
mod planetesimal_formation;
mod size_distribution;

#[cfg(test)]
mod coagulation_test;
#[cfg(test)]
mod particle_bin_test;
#[cfg(test)]
mod particle_test;
#[cfg(test)]
mod size_distribution_test;

pub use coagulation::{Coagulation, EvolutionParams};
pub use collision_outcome::{CollisionOutcome, MaterialProperties};
pub use particle::{DragRegime, Particle};
pub use particle_bin::ParticleBin;
pub use planetesimal_formation::PlanetesimalFormationEvent;
pub use size_distribution::SizeDistribution;
