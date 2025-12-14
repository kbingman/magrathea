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
    detect_ejections, detect_star_collisions,
};
pub use resolution::{merge_bodies, remove_ejections, remove_star_collisions, resolve_collisions};
