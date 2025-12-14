//! Collision resolution through momentum-conserving mergers
//!
//! When bodies collide, they merge into a single body that conserves:
//! - Total mass
//! - Total momentum
//! - Volume (assuming constant density)

use crate::body::{Body, BodyId};
use crate::collisions::CollisionEvent;
use crate::state::SystemState;
use nalgebra::Point2;
use std::collections::HashSet;

/// Merge two bodies, conserving mass and momentum
///
/// The merged body conserves:
/// - Total mass: m_new = m_a + m_b
/// - Momentum: p_new = p_a + p_b
/// - Volume: V_new = V_a + V_b (constant density assumption)
///
/// In 2D, if we assume "thickness" in the z-direction, volume scales as πr²,
/// so r_new = sqrt(r_a² + r_b²)
///
/// # Arguments
///
/// * `a` - First body
/// * `b` - Second body
/// * `new_id` - ID for the merged body
///
/// # Returns
///
/// New body representing the merger
///
/// # Examples
///
/// ```
/// use nbody::collisions::merge_bodies;
/// use nbody::body::{Body, BodyId};
/// use nalgebra::{Point2, Vector2};
///
/// let a = Body {
///     id: BodyId(0),
///     mass: 1.0,
///     radius: 0.01,
///     position: Point2::new(1.0, 0.0),
///     velocity: Vector2::new(0.0, 5.0),
/// };
///
/// let b = Body {
///     id: BodyId(1),
///     mass: 1.0,
///     radius: 0.01,
///     position: Point2::new(1.1, 0.0),
///     velocity: Vector2::new(0.0, 3.0),
/// };
///
/// let merged = merge_bodies(&a, &b, BodyId(2));
///
/// // Mass is conserved
/// assert!((merged.mass - 2.0).abs() < 1e-10);
///
/// // Momentum is conserved
/// let p_initial = a.momentum() + b.momentum();
/// let p_final = merged.momentum();
/// assert!((p_final.x - p_initial.x).abs() < 1e-10);
/// assert!((p_final.y - p_initial.y).abs() < 1e-10);
/// ```
pub fn merge_bodies(a: &Body, b: &Body, new_id: BodyId) -> Body {
    let total_mass = a.mass + b.mass;

    // Center of mass position
    let pos_coords = (a.position.coords * a.mass + b.position.coords * b.mass) / total_mass;
    let position = Point2::from(pos_coords);

    // Momentum-conserving velocity
    let velocity = (a.momentum() + b.momentum()) / total_mass;

    // Volume-conserving radius (2D: area ~ r²)
    // For 3D spheres with equal density: r = (r_a³ + r_b³)^(1/3)
    // For 2D disks: r = sqrt(r_a² + r_b²)
    let radius = (a.radius.powi(2) + b.radius.powi(2)).sqrt();

    Body {
        id: new_id,
        mass: total_mass,
        radius,
        position,
        velocity,
    }
}

/// Process all collision events and merge bodies
///
/// Handles collision cascades correctly by:
/// 1. Sorting collisions by separation (closest first)
/// 2. Tracking which bodies have been consumed
/// 3. Skipping events involving already-consumed bodies
///
/// # Arguments
///
/// * `state` - System state to modify
/// * `events` - List of collision events
///
/// # Examples
///
/// ```
/// use nbody::collisions::{resolve_collisions, CollisionEvent};
/// use nbody::body::BodyId;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
///
/// let id_a = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
/// let id_b = system.add_body(0.001, 0.01, Point2::new(1.001, 0.0), Vector2::new(0.0, 6.28));
///
/// let events = vec![CollisionEvent {
///     body_a: id_a,
///     body_b: id_b,
///     separation: 0.001,
///     collision_radius: 0.02,
/// }];
///
/// resolve_collisions(&mut system, events);
///
/// // Should have one merged body
/// assert_eq!(system.body_count(), 1);
/// ```
pub fn resolve_collisions(state: &mut SystemState, mut events: Vec<CollisionEvent>) {
    // Sort by separation (closest first) to handle cascades correctly
    events.sort_by(|a, b| a.separation.partial_cmp(&b.separation).unwrap());

    // Track which bodies have been consumed
    let mut consumed: HashSet<BodyId> = HashSet::new();

    for event in events {
        // Skip if either body was already consumed
        if consumed.contains(&event.body_a) || consumed.contains(&event.body_b) {
            continue;
        }

        // Find the bodies (they should still exist)
        let a = state.bodies.iter().find(|b| b.id == event.body_a).copied();
        let b = state.bodies.iter().find(|b| b.id == event.body_b).copied();

        if let (Some(a), Some(b)) = (a, b) {
            // Remove both bodies
            state.remove_body(event.body_a);
            state.remove_body(event.body_b);

            // Create merged body (reuse the lower ID to maintain order)
            let new_id = if event.body_a.0 < event.body_b.0 {
                event.body_a
            } else {
                event.body_b
            };
            let merged = merge_bodies(&a, &b, new_id);
            state.bodies.push(merged);

            // Mark both as consumed
            consumed.insert(event.body_a);
            consumed.insert(event.body_b);
        }
    }
}

/// Remove bodies that have collided with the star
///
/// # Arguments
///
/// * `state` - System state to modify
/// * `ids` - List of body IDs to remove
///
/// # Examples
///
/// ```
/// use nbody::collisions::resolution::remove_star_collisions;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
///
/// let id = system.add_body(0.001, 0.01, Point2::new(0.001, 0.0), Vector2::new(0.0, 50.0));
///
/// remove_star_collisions(&mut system, vec![id]);
/// assert_eq!(system.body_count(), 0);
/// ```
pub fn remove_star_collisions(state: &mut SystemState, ids: Vec<BodyId>) {
    ids.iter().for_each(|&id| {
        state.remove_body(id);
    });
}

/// Remove bodies that have escaped the system
///
/// # Arguments
///
/// * `state` - System state to modify
/// * `ids` - List of body IDs to remove
///
/// # Examples
///
/// ```
/// use nbody::collisions::resolution::remove_ejections;
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
///
/// let id = system.add_body(0.001, 0.01, Point2::new(100.0, 0.0), Vector2::new(0.0, 20.0));
///
/// remove_ejections(&mut system, vec![id]);
/// assert_eq!(system.body_count(), 0);
/// ```
pub fn remove_ejections(state: &mut SystemState, ids: Vec<BodyId>) {
    ids.iter().for_each(|&id| {
        state.remove_body(id);
    });
}
