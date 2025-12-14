//! Collision detection and resolution for N-body systems
//!
//! This module provides collision detection using various strategies
//! and resolution through momentum-conserving mergers.

pub mod detection;
pub mod resolution;

#[cfg(test)]
mod detection_test;
#[cfg(test)]
mod resolution_test;

pub use detection::{
    CollisionCriteria, CollisionDetector, CollisionEvent, DirectDetector, TreeDetector,
};
pub use resolution::{merge_bodies, resolve_collisions};
