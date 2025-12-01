//! Emergent planet formation simulation
//!
//! This module implements a bottom-up simulation where planetary system
//! architecture emerges from local physical rules rather than prescribed
//! formation pathways.

pub(crate) mod constants;
mod disk_model;
pub(crate) mod gas_disk;

#[cfg(test)]
mod disk_model_test;
#[cfg(test)]
mod gas_disk_test;

pub use disk_model::{DiskMass, DiskModel};
pub use gas_disk::GasDisk;
