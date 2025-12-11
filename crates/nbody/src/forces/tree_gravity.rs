//! Tree-based gravity using Barnes-Hut algorithm (O(N log N))

use crate::arena_bhtree::{BHTree, BoundingBox};
use crate::forces::{ForceModel, G};
use crate::state::SystemState;
use nalgebra::Vector2;

/// Barnes-Hut tree-based gravitational force computation
///
/// Uses a quadtree to approximate distant forces, achieving O(N log N)
/// complexity instead of O(N²). Accuracy controlled by the opening angle θ.
///
/// Best for:
/// - Large systems (N > 100)
/// - Production simulations
/// - Systems where performance matters
///
/// # Opening Angle (θ)
///
/// Controls the accuracy/speed tradeoff:
/// - θ = 0.0: Exact (same as DirectGravity)
/// - θ = 0.5: High accuracy, moderately fast (recommended)
/// - θ = 1.0: Lower accuracy, very fast
/// - θ = 2.0: Low accuracy, maximum speed
///
/// # Examples
///
/// ```
/// use nbody::forces::{TreeGravity, ForceModel};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
/// system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
///
/// // Use default θ = 0.5
/// let gravity = TreeGravity::new();
/// let accel = gravity.acceleration(0, &system);
/// ```
pub struct TreeGravity {
    /// Opening angle parameter (typically 0.5)
    pub theta: f64,
    /// Optional softening length (AU)
    pub softening: f64,
}

impl TreeGravity {
    /// Creates a new tree gravity force with default parameters
    ///
    /// Default θ = 0.5 (good accuracy/speed balance)
    pub fn new() -> Self {
        Self {
            theta: 0.5,
            softening: 0.0,
        }
    }

    /// Creates a tree gravity force with custom opening angle
    ///
    /// # Arguments
    ///
    /// * `theta` - Opening angle (0.0 = exact, larger = faster but less accurate)
    pub fn with_theta(theta: f64) -> Self {
        Self {
            theta,
            softening: 0.0,
        }
    }

    /// Creates a tree gravity force with softening
    ///
    /// # Arguments
    ///
    /// * `theta` - Opening angle
    /// * `softening` - Softening length in AU
    pub fn with_softening(theta: f64, softening: f64) -> Self {
        Self { theta, softening }
    }
}

impl Default for TreeGravity {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceModel for TreeGravity {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64> {
        if state.bodies.is_empty() {
            return Vector2::zeros();
        }

        let body = &state.bodies[idx];
        let eps2 = self.softening * self.softening;

        // Acceleration from star (at origin)
        let r_star = body.position.coords;
        let r2_star = r_star.magnitude_squared() + eps2;
        let r_star_mag = r2_star.sqrt();
        let star_mass = state.star.mass.to_solar_masses();
        let a_star = -r_star * (G * star_mass / (r2_star * r_star_mag));

        // Build tree for body-body forces
        let bounds = BoundingBox::new_from_bodies(&state.bodies);
        let tree = BHTree::build(&state.bodies, bounds);

        // Use tree to compute acceleration from all bodies
        // The tree includes the body itself, so we need to subtract self-force
        let a_all_bodies = tree.acceleration(body.position, self.theta);

        // Self-force would be infinite (r=0), but the tree handles this by
        // returning zero for the leaf containing this body
        let a_bodies = a_all_bodies;

        a_star + a_bodies
    }

    fn potential_energy(&self, state: &SystemState) -> f64 {
        // For potential energy, we use direct computation to ensure accuracy
        // Tree methods can introduce errors in energy conservation
        let eps2 = self.softening * self.softening;
        let star_mass = state.star.mass.to_solar_masses();

        // Star-body potential
        let star_potential: f64 = state
            .bodies
            .iter()
            .map(|b| {
                let r = (b.position.coords.magnitude_squared() + eps2).sqrt();
                -G * star_mass * b.mass / r
            })
            .sum();

        // Body-body potential (direct N²)
        let body_potential: f64 = state
            .bodies
            .iter()
            .enumerate()
            .flat_map(|(i, a)| {
                state.bodies[i + 1..].iter().map(move |b| {
                    let dr = a.position - b.position; // Vector2
                    let r = (dr.magnitude_squared() + eps2).sqrt();
                    -G * a.mass * b.mass / r
                })
            })
            .sum();

        star_potential + body_potential
    }
}
