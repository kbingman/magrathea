//! Collision detection using various strategies
//!
//! Provides multiple collision detection implementations:
//! - DirectDetector: O(N²) for small systems
//! - TreeDetector: O(N log N) using Barnes-Hut quadtree

use crate::arena_bhtree::{BHTree, BoundingBox};
use crate::body::{Body, BodyId};
use crate::forces::G;
use crate::state::SystemState;
use std::collections::HashSet;

/// Criteria for what counts as a collision
///
/// Collisions can be detected based on:
/// - Physical contact (sum of radii)
/// - Hill sphere overlap (gravitational sphere of influence)
///
/// # Examples
///
/// ```
/// use nbody::collisions::CollisionCriteria;
///
/// // Default: 0.5 Hill radii, physical collisions enabled
/// let criteria = CollisionCriteria::default();
///
/// // Only physical collisions
/// let physical_only = CollisionCriteria {
///     hill_fraction: 0.0,
///     physical_collision: true,
/// };
///
/// // Wide Hill sphere detection
/// let wide = CollisionCriteria {
///     hill_fraction: 1.0,
///     physical_collision: true,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct CollisionCriteria {
    /// Fraction of mutual Hill radius for collision
    ///
    /// Typical values:
    /// - 0.0: Only physical contact counts
    /// - 0.5: Half a Hill radius (conservative)
    /// - 1.0: Full Hill radius (aggressive)
    pub hill_fraction: f64,

    /// Whether physical contact (sum of radii) triggers collision
    pub physical_collision: bool,
}

impl Default for CollisionCriteria {
    fn default() -> Self {
        Self {
            hill_fraction: 0.5,
            physical_collision: true,
        }
    }
}

/// A detected collision between two bodies
#[derive(Debug, Clone)]
pub struct CollisionEvent {
    /// First body ID
    pub body_a: BodyId,
    /// Second body ID
    pub body_b: BodyId,
    /// Current separation distance (AU)
    pub separation: f64,
    /// Collision threshold that was exceeded (AU)
    pub collision_radius: f64,
}

/// Compute Hill radius for a body orbiting a star
///
/// The Hill radius is the approximate radius of gravitational influence
/// of a smaller body in orbit around a larger body.
///
/// # Arguments
///
/// * `mass` - Body mass in solar masses
/// * `orbital_radius` - Orbital radius in AU
/// * `star_mass` - Star mass in solar masses
///
/// # Returns
///
/// Hill radius in AU
///
/// # Examples
///
/// ```
/// use nbody::collisions::detection::hill_radius;
/// use units::Mass;
///
/// // Earth's Hill radius at 1 AU
/// let earth_mass = Mass::from_earth_masses(1.0).to_solar_masses();
/// let r_hill = hill_radius(earth_mass, 1.0, 1.0);
/// assert!((r_hill - 0.01).abs() < 0.001); // About 0.01 AU
/// ```
pub fn hill_radius(mass: f64, orbital_radius: f64, star_mass: f64) -> f64 {
    orbital_radius * (mass / (3.0 * star_mass)).powf(1.0 / 3.0)
}

/// Compute mutual Hill radius of two bodies
///
/// Uses the average of the two orbital radii and sum of masses.
///
/// # Arguments
///
/// * `body1_mass` - First body mass in solar masses
/// * `body1_radius` - First body orbital radius in AU
/// * `body2_mass` - Second body mass in solar masses
/// * `body2_radius` - Second body orbital radius in AU
/// * `star_mass` - Star mass in solar masses
///
/// # Returns
///
/// Mutual Hill radius in AU
///
/// # Examples
///
/// ```
/// use nbody::collisions::detection::mutual_hill_radius;
/// use units::Mass;
///
/// // Two Earth-mass planets at 1 and 1.5 AU
/// let m = Mass::from_earth_masses(1.0).to_solar_masses();
/// let r_hill = mutual_hill_radius(m, 1.0, m, 1.5, 1.0);
/// assert!(r_hill > 0.0);
/// ```
pub fn mutual_hill_radius(
    body1_mass: f64,
    body1_radius: f64,
    body2_mass: f64,
    body2_radius: f64,
    star_mass: f64,
) -> f64 {
    let a_avg = (body1_radius + body2_radius) / 2.0;
    let m_sum = body1_mass + body2_mass;
    a_avg * (m_sum / (3.0 * star_mass)).powf(1.0 / 3.0)
}

/// Check if a pair of bodies should collide
fn check_pair(
    a: &Body,
    b: &Body,
    star_mass: f64,
    criteria: &CollisionCriteria,
) -> Option<CollisionEvent> {
    let separation = a.distance_to(b);

    // Physical collision radius
    let physical_radius = a.radius + b.radius;

    // Hill sphere collision radius
    let r_hill = mutual_hill_radius(
        a.mass,
        a.orbital_radius(),
        b.mass,
        b.orbital_radius(),
        star_mass,
    );
    let hill_collision_radius = r_hill * criteria.hill_fraction;

    // Effective collision radius is max of physical and Hill
    let collision_radius = if criteria.physical_collision {
        physical_radius.max(hill_collision_radius)
    } else {
        hill_collision_radius
    };

    if separation < collision_radius {
        Some(CollisionEvent {
            body_a: a.id,
            body_b: b.id,
            separation,
            collision_radius,
        })
    } else {
        None
    }
}

/// Collision detector trait
///
/// Different implementations offer tradeoffs between accuracy and performance.
pub trait CollisionDetector: Send + Sync {
    /// Detect all collisions in the system
    ///
    /// # Arguments
    ///
    /// * `state` - Current system state
    /// * `criteria` - Collision criteria
    ///
    /// # Returns
    ///
    /// List of collision events
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent>;
}

/// Direct O(N²) collision detector
///
/// Checks every pair of bodies. Simple and exact, but slow for large N.
///
/// Best for:
/// - Small systems (N < 100)
/// - Testing and validation
/// - When exact collision detection is critical
///
/// # Examples
///
/// ```
/// use nbody::collisions::{CollisionDetector, DirectDetector, CollisionCriteria};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
/// system.add_body(0.001, 0.01, Point2::new(1.001, 0.0), Vector2::new(0.0, 6.28));
///
/// let detector = DirectDetector;
/// let criteria = CollisionCriteria::default();
/// let collisions = detector.detect(&system, &criteria);
/// assert!(collisions.len() > 0); // Bodies are very close
/// ```
pub struct DirectDetector;

impl CollisionDetector for DirectDetector {
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent> {
        let n = state.bodies.len();
        let star_mass = state.star.mass.to_solar_masses();

        (0..n)
            .flat_map(|i| {
                ((i + 1)..n).filter_map(move |j| {
                    check_pair(&state.bodies[i], &state.bodies[j], star_mass, criteria)
                })
            })
            .collect()
    }
}

/// Tree-based collision detector using Barnes-Hut quadtree
///
/// Uses the Barnes-Hut quadtree to efficiently find nearby bodies.
/// Most efficient when the tree is already built for gravity calculations.
///
/// Best for:
/// - Large systems (N > 100)
/// - Production simulations
/// - When tree is already available from gravity calculation
///
/// # Examples
///
/// ```
/// use nbody::collisions::{CollisionDetector, TreeDetector, CollisionCriteria};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
/// system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
///
/// let detector = TreeDetector;
/// let criteria = CollisionCriteria::default();
/// let collisions = detector.detect(&system, &criteria);
/// ```
pub struct TreeDetector;

impl CollisionDetector for TreeDetector {
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent> {
        if state.bodies.is_empty() {
            return Vec::new();
        }

        // Build tree for proximity queries
        let bounds = BoundingBox::new_from_bodies(&state.bodies);
        let tree = BHTree::build(&state.bodies, bounds);

        let star_mass = state.star.mass.to_solar_masses();
        let mut events = Vec::new();
        let mut seen: HashSet<(u32, u32)> = HashSet::new();

        for body in &state.bodies {
            // Search radius: conservative estimate of collision zone
            // Use 3x Hill radius to ensure we catch all potential collisions
            let r_hill = hill_radius(body.mass, body.orbital_radius(), star_mass);
            let search_radius = (r_hill * 3.0).max(body.radius * 10.0);

            let neighbors = tree.neighbors_within(body.position, search_radius);

            for &neighbor_idx in &neighbors {
                let neighbor = &state.bodies[neighbor_idx];

                // Skip self
                if neighbor.id == body.id {
                    continue;
                }

                // Canonical ordering to avoid duplicates
                let pair = if body.id.0 < neighbor.id.0 {
                    (body.id.0, neighbor.id.0)
                } else {
                    (neighbor.id.0, body.id.0)
                };

                if seen.contains(&pair) {
                    continue;
                }
                seen.insert(pair);

                // Check if they actually collide
                if let Some(event) = check_pair(body, neighbor, star_mass, criteria) {
                    events.push(event);
                }
            }
        }

        events
    }
}

/// Check for bodies hitting the star
///
/// Detects bodies whose orbital radius is less than the star's radius.
///
/// # Arguments
///
/// * `state` - Current system state
///
/// # Returns
///
/// List of body IDs that have collided with the star
///
/// # Examples
///
/// ```
/// use nbody::collisions::detection::detect_star_collisions;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
///
/// // Add a body very close to the star
/// system.add_body(0.001, 0.01, Point2::new(0.001, 0.0), Vector2::new(0.0, 50.0));
///
/// let star_collisions = detect_star_collisions(&system);
/// assert!(star_collisions.len() > 0);
/// ```
pub fn detect_star_collisions(state: &SystemState) -> Vec<BodyId> {
    let star_radius_au = state.star.radius.as_length().to_au();
    state
        .bodies
        .iter()
        .filter(|b| b.orbital_radius() < star_radius_au)
        .map(|b| b.id)
        .collect()
}

/// Check for bodies escaping the system
///
/// Detects bodies with positive total energy (escape velocity) beyond
/// a specified radius.
///
/// # Arguments
///
/// * `state` - Current system state
/// * `escape_radius` - Minimum radius to check for ejection (AU)
///
/// # Returns
///
/// List of body IDs that have escaped
///
/// # Examples
///
/// ```
/// use nbody::collisions::detection::detect_ejections;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
///
/// // Add a body moving at escape velocity
/// system.add_body(0.001, 0.01, Point2::new(100.0, 0.0), Vector2::new(0.0, 20.0));
///
/// let ejections = detect_ejections(&system, 50.0);
/// assert!(ejections.len() > 0);
/// ```
pub fn detect_ejections(state: &SystemState, escape_radius: f64) -> Vec<BodyId> {
    let star_mass = state.star.mass.to_solar_masses();

    state
        .bodies
        .iter()
        .filter(|b| {
            let r = b.orbital_radius();
            let v = b.orbital_velocity();
            let mu = G * star_mass;

            // Escape if: v² > 2μ/r (positive energy) AND r > escape_radius
            let escape_velocity_squared = 2.0 * mu / r;
            let v_squared = v * v;

            r > escape_radius && v_squared > escape_velocity_squared
        })
        .map(|b| b.id)
        .collect()
}
