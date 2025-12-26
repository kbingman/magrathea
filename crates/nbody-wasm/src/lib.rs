//! WASM bindings for N-body simulation.
//!
//! This crate provides JavaScript/TypeScript bindings for the nbody crate,
//! enabling browser-based animation of stellar systems.
//!
//! # Architecture
//!
//! The simulation state is stored in thread-local storage (WASM is single-threaded).
//! Functions return opaque IDs for referencing mutable state, and serializable
//! snapshots for reading state.
//!
//! # Example Usage (JavaScript)
//!
//! ```javascript
//! // Create a simulation with a solar-mass star
//! const simId = simulation_create(1.0, 0.0, 4600.0);
//!
//! // Add Earth at 1 AU with circular velocity
//! const earthId = simulation_add_body(simId, {
//!   massEarth: 1.0,
//!   radiusEarth: 1.0,
//!   position: [1.0, 0.0],
//!   velocity: [0.0, 6.28]
//! });
//!
//! // Step simulation forward
//! simulation_step(simId, 0.01);
//!
//! // Get all body positions for rendering
//! const bodies = simulation_get_bodies(simId);
//! ```

use std::cell::RefCell;
use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use nbody::body::BodyId;
use nbody::collisions::{
    CollisionCriteria, CollisionDetector, CollisionEvent, DirectDetector, TreeDetector,
};
use nbody::forces::{DirectGravity, ForceModel, TreeGravity};
use nbody::integrator::{Integrator, Leapfrog};
use nbody::state::SystemState;
use stellar::generation::main_sequence_star;
use units::{Length, Mass};

// =============================================================================
// Serialization helpers
// =============================================================================

fn to_js<T: Serialize>(value: &T) -> Result<JsValue, JsError> {
    serde_wasm_bindgen::to_value(value).map_err(|e| JsError::new(&e.to_string()))
}

fn from_js<T: serde::de::DeserializeOwned>(value: JsValue) -> Result<T, JsError> {
    serde_wasm_bindgen::from_value(value).map_err(|e| JsError::new(&e.to_string()))
}

// =============================================================================
// Thread-local storage for simulation state
// =============================================================================

/// Internal simulation state including the system and configuration
struct Simulation {
    state: SystemState,
    integrator: Leapfrog,
    force: ForceType,
    collision_detector: DetectorType,
}

enum ForceType {
    Direct(DirectGravity),
    Tree(TreeGravity),
}

impl ForceType {
    fn as_force_model(&self) -> &dyn ForceModel {
        match self {
            ForceType::Direct(g) => g,
            ForceType::Tree(g) => g,
        }
    }
}

enum DetectorType {
    Direct(DirectDetector),
    Tree(TreeDetector),
}

impl DetectorType {
    fn as_detector(&self) -> &dyn CollisionDetector {
        match self {
            DetectorType::Direct(d) => d,
            DetectorType::Tree(d) => d,
        }
    }
}

thread_local! {
    static SIMULATIONS: RefCell<HashMap<u32, Simulation>> = RefCell::new(HashMap::new());
    static NEXT_SIM_ID: RefCell<u32> = const { RefCell::new(0) };
}

// =============================================================================
// Serializable types for JavaScript interop
// =============================================================================

/// 2D vector representation for JavaScript
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

/// Body data for adding new bodies
#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BodyInput {
    /// Mass in Earth masses
    pub mass_earth: f64,
    /// Physical radius in Earth radii
    pub radius_earth: f64,
    /// Position [x, y] in AU
    pub position: [f64; 2],
    /// Velocity [x, y] in AU/year
    pub velocity: [f64; 2],
}

/// Body state snapshot for JavaScript
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BodySnapshot {
    /// Unique body ID
    pub id: u32,
    /// Mass in solar masses
    pub mass_solar: f64,
    /// Mass in Earth masses
    pub mass_earth: f64,
    /// Physical radius in AU
    pub radius_au: f64,
    /// Position [x, y] in AU
    pub position: [f64; 2],
    /// Velocity [x, y] in AU/year
    pub velocity: [f64; 2],
    /// Distance from star in AU
    pub orbital_radius: f64,
    /// Speed in AU/year
    pub orbital_velocity: f64,
}

/// Complete simulation state snapshot
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulationSnapshot {
    /// Current simulation time in years
    pub time: f64,
    /// Star mass in solar masses
    pub star_mass: f64,
    /// All bodies in the system
    pub bodies: Vec<BodySnapshot>,
    /// Total kinetic energy (code units)
    pub kinetic_energy: f64,
    /// Total potential energy (code units)
    pub potential_energy: f64,
    /// Total energy (should be conserved)
    pub total_energy: f64,
}

/// Collision event for JavaScript
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CollisionSnapshot {
    /// First body ID
    pub body_a: u32,
    /// Second body ID
    pub body_b: u32,
    /// Current separation in AU
    pub separation: f64,
    /// Collision threshold in AU
    pub collision_radius: f64,
}

impl From<&CollisionEvent> for CollisionSnapshot {
    fn from(event: &CollisionEvent) -> Self {
        Self {
            body_a: event.body_a.0,
            body_b: event.body_b.0,
            separation: event.separation,
            collision_radius: event.collision_radius,
        }
    }
}

/// Configuration for force model
#[derive(Clone, Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ForceConfig {
    /// "direct" or "tree"
    #[serde(default = "default_force_type")]
    pub force_type: String,
    /// Softening length in AU (prevents singularities)
    #[serde(default)]
    pub softening: f64,
    /// Opening angle for tree gravity (0.5 is typical)
    #[serde(default = "default_theta")]
    pub theta: f64,
}

fn default_force_type() -> String {
    "direct".to_string()
}

fn default_theta() -> f64 {
    0.5
}

impl Default for ForceConfig {
    fn default() -> Self {
        Self {
            force_type: default_force_type(),
            softening: 0.0,
            theta: default_theta(),
        }
    }
}

// =============================================================================
// Simulation management functions
// =============================================================================

/// Create a new N-body simulation.
///
/// Returns a simulation ID for use with other functions.
///
/// # Arguments
/// * `star_mass` - Star mass in solar masses
/// * `star_metallicity` - Star metallicity [Fe/H]
/// * `star_age_myr` - Star age in million years
#[wasm_bindgen]
pub fn simulation_create(star_mass: f64, star_metallicity: f64, star_age_myr: f64) -> u32 {
    let star = main_sequence_star(star_mass, star_metallicity, star_age_myr);
    let state = SystemState::new(star);

    let simulation = Simulation {
        state,
        integrator: Leapfrog::new(),
        force: ForceType::Direct(DirectGravity::new()),
        collision_detector: DetectorType::Direct(DirectDetector),
    };

    let id = NEXT_SIM_ID.with(|next_id| {
        let mut id = next_id.borrow_mut();
        let current = *id;
        *id += 1;
        current
    });

    SIMULATIONS.with(|sims| {
        sims.borrow_mut().insert(id, simulation);
    });

    id
}

/// Create a simulation with custom force configuration.
///
/// # Arguments
/// * `star_mass` - Star mass in solar masses
/// * `star_metallicity` - Star metallicity [Fe/H]
/// * `star_age_myr` - Star age in million years
/// * `config` - Force configuration object
#[wasm_bindgen]
pub fn simulation_create_with_config(
    star_mass: f64,
    star_metallicity: f64,
    star_age_myr: f64,
    config: JsValue,
) -> Result<u32, JsError> {
    let config: ForceConfig = from_js(config)?;
    let star = main_sequence_star(star_mass, star_metallicity, star_age_myr);
    let state = SystemState::new(star);

    let force = match config.force_type.as_str() {
        "tree" => ForceType::Tree(TreeGravity::with_softening(config.theta, config.softening)),
        _ => ForceType::Direct(DirectGravity::with_softening(config.softening)),
    };

    let collision_detector = match config.force_type.as_str() {
        "tree" => DetectorType::Tree(TreeDetector),
        _ => DetectorType::Direct(DirectDetector),
    };

    let simulation = Simulation {
        state,
        integrator: Leapfrog::new(),
        force,
        collision_detector,
    };

    let id = NEXT_SIM_ID.with(|next_id| {
        let mut id = next_id.borrow_mut();
        let current = *id;
        *id += 1;
        current
    });

    SIMULATIONS.with(|sims| {
        sims.borrow_mut().insert(id, simulation);
    });

    Ok(id)
}

/// Delete a simulation to free memory.
#[wasm_bindgen]
pub fn simulation_delete(sim_id: u32) {
    SIMULATIONS.with(|sims| {
        sims.borrow_mut().remove(&sim_id);
    });
}

// =============================================================================
// Body management
// =============================================================================

/// Add a body to the simulation.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `body` - Body parameters (massEarth, radiusEarth, position, velocity)
///
/// # Returns
/// The new body's ID
#[wasm_bindgen]
pub fn simulation_add_body(sim_id: u32, body: JsValue) -> Result<u32, JsError> {
    let body: BodyInput = from_js(body)?;

    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let mass = Mass::from_earth_masses(body.mass_earth).to_solar_masses();
        let radius = Length::from_earth_radii(body.radius_earth).to_au();
        let position = nalgebra::Point2::new(body.position[0], body.position[1]);
        let velocity = nalgebra::Vector2::new(body.velocity[0], body.velocity[1]);

        let id = sim.state.add_body(mass, radius, position, velocity);
        Ok(id.0)
    })
}

/// Add a body using simple parameters.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `mass_earth` - Mass in Earth masses
/// * `radius_earth` - Radius in Earth radii
/// * `x` - X position in AU
/// * `y` - Y position in AU
/// * `vx` - X velocity in AU/year
/// * `vy` - Y velocity in AU/year
///
/// # Returns
/// The new body's ID
#[wasm_bindgen]
pub fn simulation_add_body_simple(
    sim_id: u32,
    mass_earth: f64,
    radius_earth: f64,
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
) -> Result<u32, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let mass = Mass::from_earth_masses(mass_earth).to_solar_masses();
        let radius = Length::from_earth_radii(radius_earth).to_au();
        let position = nalgebra::Point2::new(x, y);
        let velocity = nalgebra::Vector2::new(vx, vy);

        let id = sim.state.add_body(mass, radius, position, velocity);
        Ok(id.0)
    })
}

/// Remove a body from the simulation.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `body_id` - Body ID to remove
///
/// # Returns
/// true if the body was found and removed
#[wasm_bindgen]
pub fn simulation_remove_body(sim_id: u32, body_id: u32) -> Result<bool, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        Ok(sim.state.remove_body(BodyId(body_id)).is_some())
    })
}

// =============================================================================
// Simulation stepping
// =============================================================================

/// Advance the simulation by one timestep.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `dt` - Timestep in years
///
/// # Returns
/// Current simulation time in years
#[wasm_bindgen]
pub fn simulation_step(sim_id: u32, dt: f64) -> Result<f64, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        sim.integrator
            .step(&mut sim.state, dt, sim.force.as_force_model());
        Ok(sim.state.time)
    })
}

/// Advance the simulation by multiple timesteps.
///
/// More efficient than calling simulation_step repeatedly.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `dt` - Timestep in years
/// * `n_steps` - Number of steps to take
///
/// # Returns
/// Current simulation time in years
#[wasm_bindgen]
pub fn simulation_integrate(sim_id: u32, dt: f64, n_steps: usize) -> Result<f64, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        sim.integrator
            .integrate(&mut sim.state, dt, n_steps, sim.force.as_force_model());
        Ok(sim.state.time)
    })
}

// =============================================================================
// State queries
// =============================================================================

/// Get all body positions and velocities for rendering.
///
/// This is optimized for animation: returns minimal data in a flat format.
///
/// # Returns
/// Array of body snapshots with positions, velocities, etc.
#[wasm_bindgen]
pub fn simulation_get_bodies(sim_id: u32) -> Result<JsValue, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let bodies: Vec<BodySnapshot> = sim
            .state
            .bodies
            .iter()
            .map(|b| BodySnapshot {
                id: b.id.0,
                mass_solar: b.mass,
                mass_earth: Mass::from_solar_masses(b.mass).to_earth_masses(),
                radius_au: b.radius,
                position: [b.position.x, b.position.y],
                velocity: [b.velocity.x, b.velocity.y],
                orbital_radius: b.orbital_radius(),
                orbital_velocity: b.orbital_velocity(),
            })
            .collect();

        to_js(&bodies)
    })
}

/// Get complete simulation state snapshot.
///
/// Includes time, bodies, and energy for diagnostics.
#[wasm_bindgen]
pub fn simulation_get_state(sim_id: u32) -> Result<JsValue, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let kinetic_energy = sim.state.kinetic_energy();
        let potential_energy = sim.force.as_force_model().potential_energy(&sim.state);

        let bodies: Vec<BodySnapshot> = sim
            .state
            .bodies
            .iter()
            .map(|b| BodySnapshot {
                id: b.id.0,
                mass_solar: b.mass,
                mass_earth: Mass::from_solar_masses(b.mass).to_earth_masses(),
                radius_au: b.radius,
                position: [b.position.x, b.position.y],
                velocity: [b.velocity.x, b.velocity.y],
                orbital_radius: b.orbital_radius(),
                orbital_velocity: b.orbital_velocity(),
            })
            .collect();

        let snapshot = SimulationSnapshot {
            time: sim.state.time,
            star_mass: sim.state.star.mass.to_solar_masses(),
            bodies,
            kinetic_energy,
            potential_energy,
            total_energy: kinetic_energy + potential_energy,
        };

        to_js(&snapshot)
    })
}

/// Get current simulation time in years.
#[wasm_bindgen]
pub fn simulation_get_time(sim_id: u32) -> Result<f64, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        Ok(sim.state.time)
    })
}

/// Get body count.
#[wasm_bindgen]
pub fn simulation_body_count(sim_id: u32) -> Result<usize, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        Ok(sim.state.body_count())
    })
}

// =============================================================================
// Collision detection and resolution
// =============================================================================

/// Detect collisions between bodies.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `hill_fraction` - Fraction of mutual Hill radius for collision (0.5 typical)
/// * `physical` - Whether to include physical contact collisions
///
/// # Returns
/// Array of collision events
#[wasm_bindgen]
pub fn simulation_detect_collisions(
    sim_id: u32,
    hill_fraction: f64,
    physical: bool,
) -> Result<JsValue, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let criteria = CollisionCriteria {
            hill_fraction,
            physical_collision: physical,
        };

        let events = sim
            .collision_detector
            .as_detector()
            .detect(&sim.state, &criteria);
        let snapshots: Vec<CollisionSnapshot> =
            events.iter().map(CollisionSnapshot::from).collect();

        to_js(&snapshots)
    })
}

/// Resolve collisions by merging bodies.
///
/// Bodies that collide are merged, conserving mass and momentum.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `hill_fraction` - Fraction of mutual Hill radius for collision
/// * `physical` - Whether to include physical contact collisions
///
/// # Returns
/// Number of mergers that occurred
#[wasm_bindgen]
pub fn simulation_resolve_collisions(
    sim_id: u32,
    hill_fraction: f64,
    physical: bool,
) -> Result<usize, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let criteria = CollisionCriteria {
            hill_fraction,
            physical_collision: physical,
        };

        let events = sim
            .collision_detector
            .as_detector()
            .detect(&sim.state, &criteria);
        let count = events.len();

        nbody::collisions::resolve_collisions(&mut sim.state, events);

        Ok(count)
    })
}

/// Detect and remove bodies that have collided with the star.
///
/// # Returns
/// Number of bodies removed
#[wasm_bindgen]
pub fn simulation_remove_star_collisions(sim_id: u32) -> Result<usize, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let ids = nbody::collisions::detect_star_collisions(&sim.state);
        let count = ids.len();

        nbody::collisions::resolution::remove_star_collisions(&mut sim.state, ids);

        Ok(count)
    })
}

/// Detect and remove bodies that have escaped the system.
///
/// # Arguments
/// * `sim_id` - Simulation ID
/// * `escape_radius` - Minimum radius to check for ejection (AU)
///
/// # Returns
/// Number of bodies removed
#[wasm_bindgen]
pub fn simulation_remove_ejections(sim_id: u32, escape_radius: f64) -> Result<usize, JsError> {
    SIMULATIONS.with(|sims| {
        let mut sims = sims.borrow_mut();
        let sim = sims
            .get_mut(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        let ids = nbody::collisions::detect_ejections(&sim.state, escape_radius);
        let count = ids.len();

        nbody::collisions::resolution::remove_ejections(&mut sim.state, ids);

        Ok(count)
    })
}

// =============================================================================
// Utility functions
// =============================================================================

/// Calculate circular orbital velocity at a given radius.
///
/// # Arguments
/// * `star_mass` - Star mass in solar masses
/// * `radius_au` - Orbital radius in AU
///
/// # Returns
/// Circular velocity in AU/year
#[wasm_bindgen]
pub fn circular_velocity(star_mass: f64, radius_au: f64) -> f64 {
    // v = sqrt(G * M / r)
    // G = 4π² in AU³ M☉⁻¹ yr⁻²
    let g = 39.478417;
    (g * star_mass / radius_au).sqrt()
}

/// Calculate orbital period at a given radius.
///
/// # Arguments
/// * `star_mass` - Star mass in solar masses
/// * `radius_au` - Orbital radius in AU
///
/// # Returns
/// Orbital period in years
#[wasm_bindgen]
pub fn orbital_period(star_mass: f64, radius_au: f64) -> f64 {
    // P = 2π * sqrt(a³ / (G * M))
    // P² = 4π² * a³ / (G * M)
    // P = sqrt(a³ / M) for G = 4π²
    (radius_au.powi(3) / star_mass).sqrt()
}

/// Suggest a timestep for stable integration.
///
/// Returns dt such that the innermost body takes ~100 steps per orbit.
///
/// # Arguments
/// * `sim_id` - Simulation ID
///
/// # Returns
/// Suggested timestep in years
#[wasm_bindgen]
pub fn simulation_suggest_timestep(sim_id: u32) -> Result<f64, JsError> {
    SIMULATIONS.with(|sims| {
        let sims = sims.borrow();
        let sim = sims
            .get(&sim_id)
            .ok_or_else(|| JsError::new(&format!("Simulation {} not found", sim_id)))?;

        if sim.state.bodies.is_empty() {
            // Default: ~100 steps per year at 1 AU
            return Ok(0.01);
        }

        // Find innermost body
        let min_radius = sim
            .state
            .bodies
            .iter()
            .map(|b| b.orbital_radius())
            .fold(f64::INFINITY, f64::min);

        let star_mass = sim.state.star.mass.to_solar_masses();
        let period = orbital_period(star_mass, min_radius);

        // ~100 steps per orbit for innermost body
        Ok(period / 100.0)
    })
}
