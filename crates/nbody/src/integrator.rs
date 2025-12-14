//! Time integration for N-body systems
//!
//! This module provides integrators for evolving N-body systems forward in time.
//! The primary integrator is the symplectic leapfrog method, which conserves
//! energy and angular momentum over long timescales.

use crate::forces::ForceModel;
use crate::state::SystemState;
use nalgebra::Vector2;

/// A time integrator for N-body systems
///
/// Integrators advance the system state forward in time by computing
/// accelerations from force models and updating positions and velocities.
pub trait Integrator: Send + Sync {
    /// Advance the system by one timestep
    ///
    /// # Arguments
    ///
    /// * `state` - Current system state (modified in place)
    /// * `dt` - Timestep in years
    /// * `force` - Force model to compute accelerations
    fn step(&self, state: &mut SystemState, dt: f64, force: &dyn ForceModel);

    /// Advance the system by multiple timesteps
    ///
    /// # Arguments
    ///
    /// * `state` - Current system state (modified in place)
    /// * `dt` - Timestep in years
    /// * `n_steps` - Number of steps to take
    /// * `force` - Force model to compute accelerations
    ///
    /// # Returns
    ///
    /// Final time after integration
    fn integrate(
        &self,
        state: &mut SystemState,
        dt: f64,
        n_steps: usize,
        force: &dyn ForceModel,
    ) -> f64 {
        for _ in 0..n_steps {
            self.step(state, dt, force);
        }
        state.time
    }
}

/// Symplectic leapfrog integrator (2nd order)
///
/// The leapfrog integrator is a symplectic method that conserves energy
/// and angular momentum to machine precision over long timescales. It's
/// the standard choice for N-body simulations.
///
/// The algorithm alternates between half-step velocity updates (kicks)
/// and full-step position updates (drifts):
///
/// 1. Kick: v(t + dt/2) = v(t) + a(t) * dt/2
/// 2. Drift: x(t + dt) = x(t) + v(t + dt/2) * dt
/// 3. Kick: v(t + dt) = v(t + dt/2) + a(t + dt) * dt/2
///
/// # Examples
///
/// ```
/// use nbody::integrator::{Integrator, Leapfrog};
/// use nbody::forces::{ForceModel, DirectGravity};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
///
/// let integrator = Leapfrog::new();
/// let force = DirectGravity::new();
/// let dt = 0.01; // 0.01 years
///
/// // Advance by one timestep
/// integrator.step(&mut system, dt, &force);
/// assert!(system.time > 0.0);
/// ```
pub struct Leapfrog {
    /// Whether to use the drift-kick-drift (DKD) variant
    ///
    /// If true, uses DKD form which is better for velocity-dependent forces.
    /// If false, uses KDK form (kick-drift-kick).
    /// For gravity-only simulations, both are equivalent.
    pub use_dkd: bool,
}

impl Leapfrog {
    /// Creates a new leapfrog integrator with default settings
    ///
    /// Uses the kick-drift-kick (KDK) form by default.
    pub fn new() -> Self {
        Self { use_dkd: false }
    }

    /// Creates a leapfrog integrator using drift-kick-drift form
    pub fn new_dkd() -> Self {
        Self { use_dkd: true }
    }

    /// Perform a kick step: update velocities by half timestep
    fn kick(&self, state: &mut SystemState, dt_half: f64, force: &dyn ForceModel) {
        // Compute all accelerations functionally
        let accelerations: Vec<Vector2<f64>> = (0..state.bodies.len())
            .map(|i| force.acceleration(i, state))
            .collect();

        // Update velocities using iterator and zip
        state
            .bodies
            .iter_mut()
            .zip(accelerations.iter())
            .for_each(|(body, accel)| {
                body.velocity += accel * dt_half;
            });
    }

    /// Perform a drift step: update positions by full timestep
    fn drift(&self, state: &mut SystemState, dt: f64) {
        // Update positions using iterator
        state.bodies.iter_mut().for_each(|body| {
            body.position += body.velocity * dt;
        });
    }
}

impl Default for Leapfrog {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for Leapfrog {
    fn step(&self, state: &mut SystemState, dt: f64, force: &dyn ForceModel) {
        if self.use_dkd {
            // Drift-Kick-Drift form
            self.drift(state, dt / 2.0);
            self.kick(state, dt, force);
            self.drift(state, dt / 2.0);
        } else {
            // Kick-Drift-Kick form (default)
            self.kick(state, dt / 2.0, force);
            self.drift(state, dt);
            self.kick(state, dt / 2.0, force);
        }

        state.time += dt;
    }
}

/// Simple Euler integrator (1st order, for testing/comparison only)
///
/// The Euler integrator is the simplest time integrator but has poor
/// accuracy and doesn't conserve energy. It's included only for testing
/// and educational purposes.
///
/// **Do not use for production simulations!** Use `Leapfrog` instead.
///
/// # Examples
///
/// ```
/// use nbody::integrator::{Integrator, Euler};
/// use nbody::forces::DirectGravity;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
///
/// let integrator = Euler;
/// let force = DirectGravity::new();
///
/// integrator.step(&mut system, 0.01, &force);
/// ```
pub struct Euler;

impl Integrator for Euler {
    fn step(&self, state: &mut SystemState, dt: f64, force: &dyn ForceModel) {
        // Compute all accelerations at current state functionally
        let accelerations: Vec<Vector2<f64>> = (0..state.bodies.len())
            .map(|i| force.acceleration(i, state))
            .collect();

        // Update positions and velocities using iterator
        state
            .bodies
            .iter_mut()
            .zip(accelerations.iter())
            .for_each(|(body, accel)| {
                body.position += body.velocity * dt;
                body.velocity += accel * dt;
            });

        state.time += dt;
    }
}
