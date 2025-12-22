//! Emergent planet formation simulation
//!
//! This module implements a bottom-up simulation where planetary system
//! architecture emerges from local physical rules rather than prescribed
//! formation pathways.

pub(crate) mod constants;
mod disk_model;
pub(crate) mod gas_disk;
mod grid_disk;
mod photoevaporation;

#[cfg(test)]
mod disk_model_test;
#[cfg(test)]
mod gas_disk_test;
#[cfg(test)]
mod grid_disk_test;
#[cfg(test)]
mod photoevaporation_test;

pub use disk_model::{DiskMass, DiskModel};
pub use gas_disk::GasDisk;
pub use grid_disk::GridDisk;
pub use photoevaporation::{PhotoevaporatingDisk, PhotoevaporationModel};
