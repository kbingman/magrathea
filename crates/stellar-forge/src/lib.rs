//! Stellar object generation for the Magrathea planetary formation system.
//!
//! This crate provides types and generation functions for stars and stellar remnants.

pub mod disk;
pub mod particles;

// Re-export stellar crate
pub use stellar;

// Re-export types for convenience
pub use stellar::{
    BlackHole, GiantStar, LuminosityClass, MainSequenceStar, NeutronStar, SpectralType,
    StellarColor, StellarObject, VariabilityType, WhiteDwarf, WhiteDwarfType,
};

// Re-export generation functions
pub use stellar::{
    black_hole, giant_star, main_sequence_star, neutron_star, sample_main_sequence_star,
    sample_stellar_object, solar_analog, stellar_object, white_dwarf,
};

pub use disk::{DiskMass, DiskModel, GasDisk, GridDisk};
pub use particles::{DragRegime, Particle, ParticleBin, SizeDistribution};
