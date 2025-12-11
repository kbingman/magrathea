//! Force models for N-body simulations
//!
//! This module provides the `ForceModel` trait and various implementations
//! for computing gravitational and other forces on bodies.

use crate::state::SystemState;
use nalgebra::Vector2;

pub mod gravity;
pub mod tree_gravity;

#[cfg(test)]
mod gravity_test;
#[cfg(test)]
mod tree_gravity_test;

pub use gravity::DirectGravity;
pub use tree_gravity::TreeGravity;

/// Gravitational constant in AU³ M☉⁻¹ year⁻²
/// G = 4π² ≈ 39.478417
pub const G: f64 = 39.478417;

/// A source of acceleration on bodies in an N-body system
///
/// Force models compute accelerations and optionally potential energy.
/// Multiple force models can be combined using `CompositeForce`.
///
/// # Examples
///
/// ```
/// use nbody::forces::{ForceModel, DirectGravity};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
///
/// let gravity = DirectGravity::new();
/// let accel = gravity.acceleration(0, &system);
/// ```
pub trait ForceModel: Send + Sync {
    /// Compute acceleration on body at index `idx` given full system state
    ///
    /// # Arguments
    ///
    /// * `idx` - Index of the body in `state.bodies`
    /// * `state` - Current system state
    ///
    /// # Returns
    ///
    /// Acceleration vector in AU/year²
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64>;

    /// Compute potential energy contribution (optional)
    ///
    /// Default implementation returns 0.0. Override for force models
    /// that contribute to potential energy (e.g., gravity).
    ///
    /// # Arguments
    ///
    /// * `state` - Current system state
    ///
    /// # Returns
    ///
    /// Potential energy in M☉ AU² year⁻²
    fn potential_energy(&self, _state: &SystemState) -> f64 {
        0.0
    }
}

/// Combine multiple force models into a single composite force
///
/// # Examples
///
/// ```
/// use nbody::forces::{CompositeForce, DirectGravity};
///
/// let composite = CompositeForce::new()
///     .with_force(DirectGravity::new());
///
/// // Later we could add gas drag, migration, etc.
/// ```
pub struct CompositeForce {
    models: Vec<Box<dyn ForceModel>>,
}

impl CompositeForce {
    /// Creates an empty composite force
    pub fn new() -> Self {
        Self { models: Vec::new() }
    }

    /// Adds a force model to the composite
    pub fn with_force<F: ForceModel + 'static>(mut self, force: F) -> Self {
        self.models.push(Box::new(force));
        self
    }
}

impl Default for CompositeForce {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceModel for CompositeForce {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64> {
        self.models
            .iter()
            .map(|f| f.acceleration(idx, state))
            .fold(Vector2::zeros(), |acc, a| acc + a)
    }

    fn potential_energy(&self, state: &SystemState) -> f64 {
        self.models.iter().map(|f| f.potential_energy(state)).sum()
    }
}
