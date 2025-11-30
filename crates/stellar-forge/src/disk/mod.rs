//! Emergent planet formation simulation
//!
//! This module implements a bottom-up simulation where planetary system
//! architecture emerges from local physical rules rather than prescribed
//! formation pathways.

pub(crate) mod constants;
pub(crate) mod gas_disk;

#[cfg(test)]
mod gas_disk_test;

pub use gas_disk::GasDisk;
