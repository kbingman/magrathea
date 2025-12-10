# N-Body Engine Development Plan (2D Planar Version)

## Overview

A 2D gravitational N-body simulation engine that integrates with the existing Magrathea statistical and emergent crates. Designed for:

1. **Stability testing** — Evolve statistically-generated systems to verify dynamical plausibility
2. **Late-stage evolution** — Giant impact phase, resonant chain formation
3. **What-if experiments** — Perturb systems and observe outcomes
4. **Gas disk effects** — Migration and drag in protoplanetary disks

**Architecture Decision: 2D Planar Simulation**

This implementation uses 2D coordinates (Point2, Vector2) rather than full 3D. Rationale:
- Protoplanetary disks are inherently thin (h/r ~ 0.05)
- Most planetary systems have low mutual inclinations (<5°)
- Simpler quadtree (4 children) vs octree (8 children) for Barnes-Hut
- Faster performance, better cache locality
- Good enough for disk migration, stability testing, and collision dynamics
- Can be extended to 3D later if needed for inclination dynamics

**Timeline:** 6-8 weeks for full implementation with gas effects
**Minimum viable:** 2-3 weeks for pure gravity + collisions

---

## Crate Structure

```
magrathea/
├── crates/
│   ├── stellar/           # Existing
│   ├── planetary/         # Existing
│   ├── statistical/       # Existing
│   ├── emergent/          # In progress
│   └── nbody/             # NEW
│       ├── Cargo.toml
│       └── src/
│           ├── lib.rs
│           ├── body.rs           # Body, BodyId (2D)
│           ├── state.rs          # SystemState, snapshot/restore
│           ├── arena_bhtree.rs   # Barnes-Hut quadtree (DONE)
│           ├── forces/
│           │   ├── mod.rs        # ForceModel trait, CompositeForce
│           │   ├── gravity.rs    # Direct N-body gravity
│           │   ├── tree.rs       # Tree-based gravity using arena_bhtree
│           │   └── gas/
│           │       ├── mod.rs
│           │       ├── disk.rs   # GasDiskModel trait + implementations
│           │       ├── drag.rs   # Aerodynamic drag
│           │       └── migration.rs  # Type I/II migration
│           ├── integrators/
│           │   ├── mod.rs
│           │   ├── leapfrog.rs   # Symplectic integrator
│           │   └── adaptive.rs   # Adaptive timestep wrapper
│           ├── collisions/
│           │   ├── mod.rs        # CollisionDetector trait
│           │   ├── detection.rs  # DirectDetector, RadialBinDetector, TreeDetector
│           │   └── resolution.rs # Merger physics, momentum conservation
│           ├── io/
│           │   ├── mod.rs
│           │   ├── import.rs     # From System JSON
│           │   └── export.rs     # To System JSON
│           └── diagnostics/
│               ├── mod.rs
│               ├── energy.rs     # Energy conservation tracking
│               └── stability.rs  # AMD, Hill stability metrics
```

---

## Unit System

**All internal calculations use astronomical units:**
- **Length**: AU (Astronomical Units)
- **Mass**: Solar masses (M☉)
- **Time**: Years
- **Velocity**: AU/year
- **Acceleration**: AU/year²

**Gravitational constant**: G = 4π² ≈ 39.478 AU³ M☉⁻¹ year⁻²

**Rationale:**
- Natural for planetary systems (Earth at 1 AU, Sun at 1 M☉)
- Avoids precision issues with huge numbers (e.g., 1.496×10¹¹ meters)
- Matches existing units crate design
- Easy conversion via units::Length, units::Mass, units::Velocity

---

## Phase 1: Core Infrastructure (Week 1) ✅ PARTIALLY DONE

### Goals
- Basic data structures
- Vector math via nalgebra
- Body representation with 2D state

### 1.1 Vector Type ✅ DONE

Using `nalgebra::Point2<f64>` for position and `nalgebra::Vector2<f64>` for velocity.
No need for custom Vec3 implementation.

### 1.2 Body Representation ✅ DONE

```rust
// src/body.rs

use nalgebra::{Point2, Vector2};
use units::{Length, Mass};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyId(pub u32);

#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub id: BodyId,
    pub mass: f64,              // Solar masses
    pub radius: f64,            // AU (physical radius for collisions)
    pub position: Point2<f64>,  // AU (heliocentric Cartesian, 2D)
    pub velocity: Vector2<f64>, // AU/year
}

impl Body {
    pub fn momentum(&self) -> Vector2<f64> {
        self.velocity * self.mass
    }

    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.magnitude_squared()
    }

    pub fn distance_to(&self, other: &Body) -> f64 {
        (self.position - other.position).magnitude()
    }

    pub fn orbital_radius(&self) -> f64 {
        self.position.coords.magnitude()
    }

    pub fn orbital_velocity(&self) -> f64 {
        self.velocity.magnitude()
    }

    /// Angular momentum scalar (r × v in 2D, returns z-component)
    pub fn specific_angular_momentum(&self) -> f64 {
        self.position.x * self.velocity.y - self.position.y * self.velocity.x
    }
}

// Implement Massive trait for Barnes-Hut tree
impl crate::arena_bhtree::Massive for Body {
    fn position(&self) -> Point2<f64> {
        self.position
    }

    fn mass(&self) -> Mass {
        Mass::from_solar_masses(self.mass)
    }

    fn mass_solar(&self) -> f64 {
        self.mass
    }
}
```

### 1.3 System State

```rust
// src/state.rs

use crate::body::{Body, BodyId};

#[derive(Debug, Clone)]
pub struct Star {
    pub mass: f64,      // Solar masses
    pub radius: f64,    // AU
}

#[derive(Debug, Clone)]
pub struct SystemState {
    pub time: f64,              // Years since epoch
    pub star: Star,
    pub bodies: Vec<Body>,
    next_id: u32,
}

impl SystemState {
    pub fn new(star: Star) -> Self {
        Self {
            time: 0.0,
            star,
            bodies: Vec::new(),
            next_id: 0,
        }
    }

    pub fn add_body(&mut self, mass: f64, radius: f64, position: Point2<f64>, velocity: Vector2<f64>) -> BodyId {
        let id = BodyId(self.next_id);
        self.next_id += 1;
        self.bodies.push(Body { id, mass, radius, position, velocity });
        id
    }

    pub fn remove_body(&mut self, id: BodyId) -> Option<Body> {
        self.bodies.iter()
            .position(|b| b.id == id)
            .map(|idx| self.bodies.remove(idx))
    }

    pub fn get_body(&self, id: BodyId) -> Option<&Body> {
        self.bodies.iter().find(|b| b.id == id)
    }

    pub fn get_body_mut(&mut self, id: BodyId) -> Option<&mut Body> {
        self.bodies.iter_mut().find(|b| b.id == id)
    }

    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
}
```

### 1.4 Orbital Elements Conversion (2D Version)

For 2D circular/elliptical orbits, we only need:
- Semi-major axis (a)
- Eccentricity (e)
- Argument of periapsis (ω)
- True anomaly (ν)

No inclination, no longitude of ascending node (these are 3D concepts).

```rust
// src/body.rs (continued)

use std::f64::consts::PI;

pub struct OrbitalElements2D {
    pub semi_major_axis: f64,       // AU
    pub eccentricity: f64,
    pub argument_of_periapsis: f64, // radians
    pub true_anomaly: f64,          // radians
}

impl Body {
    pub fn orbital_elements(&self, mu: f64) -> OrbitalElements2D {
        let r = self.position.coords;
        let v = self.velocity;
        let r_mag = r.magnitude();
        let v_mag = v.magnitude();

        // Specific angular momentum (scalar in 2D)
        let h = self.specific_angular_momentum();

        // Eccentricity vector magnitude
        // e_vec = (v × h)/μ - r/|r|
        // In 2D: e_x = (v_y * h_z)/μ - r_x/|r|
        //        e_y = (-v_x * h_z)/μ - r_y/|r|
        let e_x = (v.y * h) / mu - r.x / r_mag;
        let e_y = (-v.x * h) / mu - r.y / r_mag;
        let e = (e_x * e_x + e_y * e_y).sqrt();

        // Semi-major axis (vis-viva)
        let energy = v_mag * v_mag / 2.0 - mu / r_mag;
        let a = -mu / (2.0 * energy);

        // Argument of periapsis (angle of eccentricity vector)
        let omega = e_y.atan2(e_x);

        // True anomaly (angle from periapsis to current position)
        let r_angle = r.y.atan2(r.x);
        let nu = (r_angle - omega + 2.0 * PI) % (2.0 * PI);

        OrbitalElements2D {
            semi_major_axis: a,
            eccentricity: e,
            argument_of_periapsis: omega,
            true_anomaly: nu,
        }
    }
}
```

### Phase 1 Deliverables

- [x] `Body` with 2D position/velocity using nalgebra types
- [x] `Body` implements `Massive` trait for Barnes-Hut tree
- [x] Constructor methods using units crate
- [ ] `SystemState` with add/remove/lookup
- [ ] `OrbitalElements2D` conversion (Cartesian ↔ Keplerian)
- [ ] Unit tests for orbital element conversion
- [ ] Test: circular orbit elements should give e ≈ 0

### Phase 1 Test

```rust
#[test]
fn circular_orbit_elements_2d() {
    const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²
    let mu = G * 1.0; // Sun mass
    let a = 1.0;  // 1 AU
    let v_circular = (mu / a).sqrt(); // Circular velocity

    let body = Body {
        id: BodyId(0),
        mass: Mass::from_earth_masses(1.0).to_solar_masses(),
        radius: Length::from_earth_radii(1.0).to_au(),
        position: Point2::new(a, 0.0),
        velocity: Vector2::new(0.0, v_circular),
    };

    let elements = body.orbital_elements(mu);

    assert!((elements.semi_major_axis - a).abs() / a < 1e-10);
    assert!(elements.eccentricity < 1e-10);
}
```

---

## Phase 2: Gravity and Integration (Week 2)

### Goals
- Force calculation
- Leapfrog integrator
- Energy conservation verification

### 2.1 Force Trait

```rust
// src/forces/mod.rs

use crate::body::Body;
use crate::state::SystemState;
use nalgebra::Vector2;

/// A source of acceleration on bodies
pub trait ForceModel: Send + Sync {
    /// Compute acceleration on body at index `idx` given full system state
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64>;

    /// Optional: compute potential energy contribution
    fn potential_energy(&self, _state: &SystemState) -> f64 {
        0.0
    }
}

/// Combine multiple force models
pub struct CompositeForce {
    models: Vec<Box<dyn ForceModel>>,
}

impl CompositeForce {
    pub fn new() -> Self {
        Self { models: Vec::new() }
    }

    pub fn add<F: ForceModel + 'static>(mut self, force: F) -> Self {
        self.models.push(Box::new(force));
        self
    }
}

impl ForceModel for CompositeForce {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64> {
        self.models.iter()
            .map(|f| f.acceleration(idx, state))
            .fold(Vector2::zeros(), |acc, a| acc + a)
    }

    fn potential_energy(&self, state: &SystemState) -> f64 {
        self.models.iter()
            .map(|f| f.potential_energy(state))
            .sum()
    }
}
```

### 2.2 Direct Gravity (2D)

```rust
// src/forces/gravity.rs

use crate::forces::ForceModel;
use crate::state::SystemState;
use nalgebra::Vector2;

pub const G: f64 = 39.478417;  // AU³ M☉⁻¹ year⁻²

pub struct DirectGravity {
    /// Optional softening length to prevent singularities
    pub softening: f64,
}

impl DirectGravity {
    pub fn new() -> Self {
        Self { softening: 0.0 }
    }

    pub fn with_softening(softening: f64) -> Self {
        Self { softening }
    }
}

impl ForceModel for DirectGravity {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64> {
        let body = &state.bodies[idx];
        let eps2 = self.softening * self.softening;

        // Acceleration from star
        let r_star = body.position.coords;  // Star at origin
        let r2_star = r_star.magnitude_squared() + eps2;
        let r_star_mag = r2_star.sqrt();
        let a_star = -r_star * (G * state.star.mass / (r2_star * r_star_mag));

        // Acceleration from other bodies
        let a_bodies = state.bodies.iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .map(|(_, other)| {
                let dr = (other.position - body.position).coords;
                let r2 = dr.magnitude_squared() + eps2;
                let r = r2.sqrt();
                dr * (G * other.mass / (r2 * r))
            })
            .fold(Vector2::zeros(), |acc, a| acc + a);

        a_star + a_bodies
    }

    fn potential_energy(&self, state: &SystemState) -> f64 {
        let eps2 = self.softening * self.softening;

        // Star-body potential
        let star_potential: f64 = state.bodies.iter()
            .map(|b| {
                let r = (b.position.coords.magnitude_squared() + eps2).sqrt();
                -G * state.star.mass * b.mass / r
            })
            .sum();

        // Body-body potential (each pair counted once)
        let body_potential: f64 = state.bodies.iter()
            .enumerate()
            .flat_map(|(i, a)| {
                state.bodies[i+1..].iter().map(move |b| {
                    let r = ((a.position - b.position).coords.magnitude_squared() + eps2).sqrt();
                    -G * a.mass * b.mass / r
                })
            })
            .sum();

        star_potential + body_potential
    }
}
```

### 2.3 Leapfrog Integrator

```rust
// src/integrators/leapfrog.rs

use crate::forces::ForceModel;
use crate::state::SystemState;
use nalgebra::Vector2;

pub struct Leapfrog<F: ForceModel> {
    force: F,
    dt: f64, // Years
}

impl<F: ForceModel> Leapfrog<F> {
    pub fn new(force: F, dt: f64) -> Self {
        Self { force, dt }
    }

    /// Single integration step (Kick-Drift-Kick)
    pub fn step(&self, state: &mut SystemState) {
        let dt = self.dt;
        let n = state.bodies.len();

        // Compute accelerations at current positions
        let accelerations: Vec<Vector2<f64>> = (0..n)
            .map(|i| self.force.acceleration(i, state))
            .collect();

        // Kick: v += a * dt/2
        for (body, &acc) in state.bodies.iter_mut().zip(&accelerations) {
            body.velocity = body.velocity + acc * (dt / 2.0);
        }

        // Drift: r += v * dt
        for body in &mut state.bodies {
            body.position += body.velocity * dt;
        }

        // Recompute accelerations at new positions
        let accelerations: Vec<Vector2<f64>> = (0..n)
            .map(|i| self.force.acceleration(i, state))
            .collect();

        // Kick: v += a * dt/2
        for (body, &acc) in state.bodies.iter_mut().zip(&accelerations) {
            body.velocity = body.velocity + acc * (dt / 2.0);
        }

        state.time += dt;
    }

    /// Run for specified duration
    pub fn evolve(&self, state: &mut SystemState, duration: f64) {
        let steps = (duration / self.dt).ceil() as usize;
        for _ in 0..steps {
            self.step(state);
        }
    }
}
```

### 2.4 Energy Diagnostics

```rust
// src/diagnostics/energy.rs

use crate::forces::ForceModel;
use crate::state::SystemState;

pub fn total_kinetic_energy(state: &SystemState) -> f64 {
    state.bodies.iter()
        .map(|b| b.kinetic_energy())
        .sum()
}

pub fn total_energy<F: ForceModel>(state: &SystemState, force: &F) -> f64 {
    total_kinetic_energy(state) + force.potential_energy(state)
}

pub struct EnergyTracker {
    pub initial_energy: f64,
    pub samples: Vec<(f64, f64)>,  // (time, relative_error)
}

impl EnergyTracker {
    pub fn new<F: ForceModel>(state: &SystemState, force: &F) -> Self {
        Self {
            initial_energy: total_energy(state, force),
            samples: Vec::new(),
        }
    }

    pub fn record<F: ForceModel>(&mut self, state: &SystemState, force: &F) {
        let current = total_energy(state, force);
        let relative_error = (current - self.initial_energy) / self.initial_energy.abs();
        self.samples.push((state.time, relative_error));
    }

    pub fn max_error(&self) -> f64 {
        self.samples.iter()
            .map(|(_, e)| e.abs())
            .fold(0.0, f64::max)
    }
}
```

### Phase 2 Deliverables

- [ ] `ForceModel` trait with `acceleration()` and `potential_energy()`
- [ ] `DirectGravity` with optional softening
- [ ] `Leapfrog` integrator
- [ ] `EnergyTracker` for conservation testing
- [ ] Test: Two-body problem traces ellipse
- [ ] Test: 10 Myr Solar System integration, energy error < 10⁻⁶

### Phase 2 Test: Two-Body Problem

```rust
#[test]
fn two_body_ellipse_2d() {
    const G: f64 = 39.478417;
    let mu = G * 1.0; // Sun
    let a = 1.0;      // 1 AU
    let e = 0.3;

    // Start at perihelion
    let r_peri = a * (1.0 - e);
    let v_peri = ((1.0 + e) * mu / (a * (1.0 - e))).sqrt();

    let star = Star { mass: 1.0, radius: Length::from_solar_radii(1.0).to_au() };
    let mut state = SystemState::new(star);
    state.add_body(
        Mass::from_earth_masses(1.0).to_solar_masses(),
        Length::from_earth_radii(1.0).to_au(),
        Point2::new(r_peri, 0.0),
        Vector2::new(0.0, v_peri),
    );

    let force = DirectGravity::new();
    let period = 2.0 * PI * (a.powi(3) / mu).sqrt();
    let dt = period / 1000.0;
    let integrator = Leapfrog::new(force, dt);

    let mut tracker = EnergyTracker::new(&state, &integrator.force);

    // Evolve for 10 orbits
    for _ in 0..10 {
        integrator.evolve(&mut state, period);
        tracker.record(&state, &integrator.force);
    }

    // Energy should be conserved to high precision
    assert!(tracker.max_error() < 1e-8);

    // Should be back near starting position
    let final_r = state.bodies[0].position.coords.magnitude();
    assert!((final_r - r_peri).abs() / r_peri < 1e-4);
}
```

---

## Phase 3: Collision Handling (Week 3)

### Goals
- Detect close approaches
- Resolve collisions via merging
- Handle edge cases (ejections, star collisions)

### 3.1 Collision Detection (2D)

The 2D version uses the existing quadtree from `arena_bhtree.rs`.

```rust
// src/collisions/detection.rs

use crate::body::{Body, BodyId};
use crate::state::SystemState;
use nalgebra::Point2;

pub const G: f64 = 39.478417;

/// Criteria for what counts as a collision
#[derive(Debug, Clone)]
pub struct CollisionCriteria {
    /// Fraction of mutual Hill radius for collision
    pub hill_fraction: f64,
    /// Minimum distance (physical radii always trigger)
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

#[derive(Debug, Clone)]
pub struct CollisionEvent {
    pub body_a: BodyId,
    pub body_b: BodyId,
    pub separation: f64,
    pub collision_radius: f64,
}

/// Mutual Hill radius of two bodies orbiting a star
pub fn mutual_hill_radius(m1: f64, m2: f64, a1: f64, a2: f64, m_star: f64) -> f64 {
    let a_avg = (a1 + a2) / 2.0;
    let m_sum = m1 + m2;
    a_avg * (m_sum / (3.0 * m_star)).powf(1.0 / 3.0)
}

pub fn hill_radius(mass: f64, orbital_radius: f64, star_mass: f64) -> f64 {
    orbital_radius * (mass / (3.0 * star_mass)).powf(1.0 / 3.0)
}

// ============================================================================
// Collision Detector Trait
// ============================================================================

pub trait CollisionDetector: Send + Sync {
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent>;
}

// ============================================================================
// Direct O(N²) - Baseline for Small N
// ============================================================================

pub struct DirectDetector;

impl CollisionDetector for DirectDetector {
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent> {
        let mut events = Vec::new();
        let n = state.bodies.len();

        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(event) = check_pair(&state.bodies[i], &state.bodies[j], state.star.mass, criteria) {
                    events.push(event);
                }
            }
        }

        events
    }
}

fn check_pair(a: &Body, b: &Body, star_mass: f64, criteria: &CollisionCriteria) -> Option<CollisionEvent> {
    let separation = a.distance_to(b);

    // Physical collision radius
    let physical_radius = a.radius + b.radius;

    // Hill sphere collision radius
    let r_hill = mutual_hill_radius(
        a.mass, b.mass,
        a.orbital_radius(), b.orbital_radius(),
        star_mass,
    );
    let hill_collision_radius = r_hill * criteria.hill_fraction;

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

// ============================================================================
// Tree-Based - Uses existing arena_bhtree
// ============================================================================

use crate::arena_bhtree::BHTree;

/// Collision detection using Barnes-Hut quadtree neighbor queries.
/// Most efficient when tree is already built for gravity calculation.
pub struct TreeDetector;

impl TreeDetector {
    /// Detect collisions using a pre-built quadtree
    pub fn detect_with_tree(
        &self,
        tree: &BHTree<Body>,
        state: &SystemState,
        criteria: &CollisionCriteria,
    ) -> Vec<CollisionEvent> {
        let mut events = Vec::new();
        let mut seen: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();

        for body in &state.bodies {
            let r_hill = hill_radius(body.mass, body.orbital_radius(), state.star.mass);
            // Search for neighbors within ~3x Hill radius
            let search_radius = r_hill * 3.0;

            let neighbors = tree.neighbors_within(body.position, search_radius);

            for neighbor_idx in neighbors {
                let neighbor = &state.bodies[neighbor_idx];
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

                if let Some(event) = check_pair(body, neighbor, state.star.mass, criteria) {
                    events.push(event);
                }
            }
        }

        events
    }
}

// ============================================================================
// Star Collision and Ejection Detection
// ============================================================================

/// Check for bodies hitting the star
pub fn detect_star_collisions(state: &SystemState) -> Vec<BodyId> {
    state.bodies.iter()
        .filter(|b| b.orbital_radius() < state.star.radius)
        .map(|b| b.id)
        .collect()
}

/// Check for bodies escaping the system
pub fn detect_ejections(state: &SystemState, escape_radius: f64) -> Vec<BodyId> {
    state.bodies.iter()
        .filter(|b| {
            let r = b.orbital_radius();
            let v = b.orbital_velocity();
            let mu = G * state.star.mass;

            // Escape if: v² > 2μ/r (positive energy) AND r > escape_radius
            let escape_velocity = (2.0 * mu / r).sqrt();
            r > escape_radius && v > escape_velocity
        })
        .map(|b| b.id)
        .collect()
}
```

### 3.2 Collision Resolution (2D)

```rust
// src/collisions/resolution.rs

use crate::body::{Body, BodyId};
use crate::state::SystemState;
use nalgebra::{Point2, Vector2};

/// Merge two bodies, conserving momentum
pub fn merge_bodies(a: &Body, b: &Body, new_id: BodyId) -> Body {
    let total_mass = a.mass + b.mass;

    // Center of mass position (convert Point2 to coords, then back)
    let pos_coords = (a.position.coords * a.mass + b.position.coords * b.mass) / total_mass;
    let position = Point2::from(pos_coords);

    // Momentum-conserving velocity
    let velocity = (a.momentum() + b.momentum()) / total_mass;

    // Combined radius (volume-conserving, assuming same density)
    // In 2D: A = πr², so r_new = sqrt(r_a² + r_b²)
    let radius = (a.radius.powi(2) + b.radius.powi(2)).sqrt();

    Body {
        id: new_id,
        mass: total_mass,
        radius,
        position,
        velocity,
    }
}

/// Process all collisions in a state, returning the modified state
pub fn resolve_collisions(state: &mut SystemState, events: Vec<super::detection::CollisionEvent>) {
    // Sort by separation (closest first) to handle cascades correctly
    let mut events = events;
    events.sort_by(|a, b| a.separation.partial_cmp(&b.separation).unwrap());

    // Track which bodies have been consumed
    let mut consumed: std::collections::HashSet<BodyId> = std::collections::HashSet::new();

    for event in events {
        // Skip if either body was already consumed
        if consumed.contains(&event.body_a) || consumed.contains(&event.body_b) {
            continue;
        }

        // Find the bodies
        let a = state.bodies.iter().find(|b| b.id == event.body_a).copied();
        let b = state.bodies.iter().find(|b| b.id == event.body_b).copied();

        if let (Some(a), Some(b)) = (a, b) {
            // Remove both bodies
            state.remove_body(event.body_a);
            state.remove_body(event.body_b);

            // Create merged body (reuse the lower ID)
            let new_id = if event.body_a.0 < event.body_b.0 { event.body_a } else { event.body_b };
            let merged = merge_bodies(&a, &b, new_id);
            state.bodies.push(merged);

            // Mark both as consumed
            consumed.insert(event.body_a);
            consumed.insert(event.body_b);
        }
    }
}

/// Remove bodies that hit the star
pub fn remove_star_collisions(state: &mut SystemState, ids: Vec<BodyId>) {
    for id in ids {
        state.remove_body(id);
    }
}

/// Remove ejected bodies
pub fn remove_ejections(state: &mut SystemState, ids: Vec<BodyId>) {
    for id in ids {
        state.remove_body(id);
    }
}
```

### Phase 3 Deliverables

- [ ] `CollisionDetector` trait with two implementations:
  - [ ] `DirectDetector` — O(N²) baseline for small N
  - [ ] `TreeDetector` — uses existing arena_bhtree
- [ ] `CollisionCriteria` configuration
- [ ] `merge_bodies()` with momentum conservation
- [ ] `detect_star_collisions()` and `detect_ejections()`
- [ ] `evolve_with_collisions()` integrated loop
- [ ] Test: Two bodies on collision course merge correctly
- [ ] Test: Body on escape trajectory is removed

---

## Phase 4: System I/O (Week 4)

### Goals
- Import from existing `System` JSON format
- Export evolved state to `System` JSON
- Round-trip validation

### 4.1 Import from System (2D Version)

```rust
// src/io/import.rs

use crate::body::Body;
use crate::state::{Star, SystemState};
use nalgebra::{Point2, Vector2};
use std::f64::consts::PI;

pub const G: f64 = 39.478417;

/// Convert 2D orbital elements to Cartesian state
pub fn elements_to_cartesian_2d(
    a: f64,           // semi-major axis (AU)
    e: f64,           // eccentricity
    omega: f64,       // argument of periapsis (rad)
    nu: f64,          // true anomaly (rad)
    mu: f64,          // G * M_star
) -> (Point2<f64>, Vector2<f64>) {
    // Distance from focus
    let r = a * (1.0 - e * e) / (1.0 + e * nu.cos());

    // Position in orbital frame
    let x_orb = r * nu.cos();
    let y_orb = r * nu.sin();

    // Velocity in orbital frame
    let h = (mu * a * (1.0 - e * e)).sqrt();  // specific angular momentum
    let vx_orb = -mu / h * nu.sin();
    let vy_orb = mu / h * (e + nu.cos());

    // Rotate by argument of periapsis
    let cos_w = omega.cos();
    let sin_w = omega.sin();

    let x = cos_w * x_orb - sin_w * y_orb;
    let y = sin_w * x_orb + cos_w * y_orb;

    let vx = cos_w * vx_orb - sin_w * vy_orb;
    let vy = sin_w * vx_orb + cos_w * vy_orb;

    (Point2::new(x, y), Vector2::new(vx, vy))
}

/// Import configuration
pub struct ImportConfig {
    /// Initial true anomaly distribution
    pub anomaly_mode: AnomalyMode,
    /// RNG seed for randomized anomalies
    pub seed: u64,
}

pub enum AnomalyMode {
    /// All bodies at periapsis
    Periapsis,
    /// All bodies at apoapsis
    Apoapsis,
    /// Random true anomalies (seeded)
    Random,
    /// Evenly distributed around orbits
    Uniform,
}

impl Default for ImportConfig {
    fn default() -> Self {
        Self {
            anomaly_mode: AnomalyMode::Random,
            seed: 42,
        }
    }
}

/// Placeholder for actual System type - replace with real import
pub struct SystemJson {
    pub star_mass: f64,       // Solar masses
    pub star_radius: f64,     // Solar radii
    pub planets: Vec<PlanetJson>,
}

pub struct PlanetJson {
    pub mass: f64,            // Earth masses
    pub radius: f64,          // Earth radii
    pub semi_major_axis: f64, // AU
    pub eccentricity: f64,
}

pub fn import_system(json: &SystemJson, config: &ImportConfig) -> SystemState {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use units::{Length, Mass};

    let star = Star {
        mass: json.star_mass,
        radius: Length::from_solar_radii(json.star_radius).to_au(),
    };

    let mut state = SystemState::new(star);
    let mut rng = ChaCha8Rng::seed_from_u64(config.seed);

    for (i, planet) in json.planets.iter().enumerate() {
        let a = planet.semi_major_axis;
        let e = planet.eccentricity;
        let mass = Mass::from_earth_masses(planet.mass).to_solar_masses();
        let radius = Length::from_earth_radii(planet.radius).to_au();

        // Determine true anomaly
        let nu = match config.anomaly_mode {
            AnomalyMode::Periapsis => 0.0,
            AnomalyMode::Apoapsis => PI,
            AnomalyMode::Random => rng.gen::<f64>() * 2.0 * PI,
            AnomalyMode::Uniform => (i as f64 / json.planets.len() as f64) * 2.0 * PI,
        };

        // Random argument of periapsis
        let omega = rng.gen::<f64>() * 2.0 * PI;

        let mu = G * star.mass;
        let (position, velocity) = elements_to_cartesian_2d(a, e, omega, nu, mu);

        state.add_body(mass, radius, position, velocity);
    }

    state
}
```

### 4.2 Export to System (2D)

```rust
// src/io/export.rs

use crate::state::SystemState;
use units::{Length, Mass};

use super::import::{SystemJson, PlanetJson};

pub fn export_system(state: &SystemState) -> SystemJson {
    use crate::forces::gravity::G;

    let mu = G * state.star.mass;

    let planets = state.bodies.iter()
        .map(|body| {
            let elements = body.orbital_elements(mu);

            PlanetJson {
                mass: Mass::from_solar_masses(body.mass).to_earth_masses(),
                radius: Length::from_au(body.radius).to_earth_radii(),
                semi_major_axis: elements.semi_major_axis,
                eccentricity: elements.eccentricity,
            }
        })
        .collect();

    SystemJson {
        star_mass: state.star.mass,
        star_radius: Length::from_au(state.star.radius).to_solar_radii(),
        planets,
    }
}
```

### Phase 4 Deliverables

- [ ] `elements_to_cartesian_2d()` conversion
- [ ] `import_system()` with configurable anomaly modes
- [ ] `export_system()` to JSON
- [ ] Round-trip test: import → export preserves orbital elements
- [ ] Integration test: import statistical output, evolve, export

---

## Phase 5: Gas Disk Effects (Week 5-6)

*Same as original plan, but using 2D coordinates. Gas disk models naturally assume thin disks anyway.*

---

## Phase 6: Barnes-Hut Quadtree ✅ DONE

The quadtree implementation in `arena_bhtree.rs` is already complete and working with 2D coordinates.

### Phase 6 Status

- [x] `BoundingBox` with quadrant subdivision
- [x] `Node` enum with Empty/Leaf/LeafMulti/Internal
- [x] Tree construction from body list
- [x] Tree-walk acceleration with opening angle θ
- [x] `neighbors_within()` for proximity queries
- [x] All tests passing

---

## Testing Strategy

### Unit Tests (Per Module)

Each module should have tests for:
- Edge cases (zero mass, zero distance, etc.)
- Known analytic solutions
- Conservation laws

### Integration Tests

1. **Two-body problem**: Verify elliptical orbit, period matches Kepler
2. **Solar System stability**: 10 Myr integration, energy conservation
3. **Collision cascade**: Many small bodies → few large bodies
4. **Migration track**: Planet migrates inward, stops when disk disperses
5. **Round-trip I/O**: Statistical → N-body → export preserves elements

### Validation Tests

Compare against:
- REBOUND (established N-body code, has 2D mode)
- Known Solar System ephemerides (projected to 2D)
- Published migration timescales

---

## Success Criteria

### Minimum Viable (Phases 1-4)

- [ ] Evolve a 10-planet system for 1 Gyr
- [ ] Energy conservation < 10⁻⁵ relative error
- [ ] Correctly identify unstable systems (planets ejected)
- [ ] Round-trip with statistical generator

### Full Implementation (Phases 1-6)

- [ ] All minimum viable criteria
- [ ] Gas disk effects produce inward migration
- [ ] Migration stops when disk disperses
- [ ] Can simulate 10,000 planetesimals (with quadtree)
- [ ] Reproduce known migration timescales within factor of 2

---

## Dependencies

### Required
- `nalgebra` = "0.34.1" - for Point2, Vector2 math
- `units` (internal) - for Mass, Length conversions
- `rand` + `rand_chacha` - for seeded RNG

### Optional
- `rayon` - for parallel force calculation
- `serde` - for JSON I/O (likely already have)

### Development
- `criterion` - for benchmarking
- `proptest` - for property-based testing

---

## Future Extensions

### If 3D becomes necessary:

1. Change `Point2` → `Point3`, `Vector2` → `Vector3` in Body
2. Change quadtree → octree in arena_bhtree
3. Add full orbital elements (i, Ω) conversion
4. Add vertical structure to gas disk model

This should be straightforward because:
- Force calculation math is dimension-agnostic
- Integrator doesn't care about dimensionality
- Collision detection algorithm is the same

### Other possible extensions:

- Adaptive timestep per body (individual timesteps)
- Relativistic corrections for close-in planets
- Tidal effects (spin-orbit coupling)
- Radiation pressure for small particles
- Binary star support (two massive bodies at origin)

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Integration instability | Use symplectic integrator, verify energy conservation |
| Close encounter singularities | Collision detection + merging before singularity |
| Gas physics too complex | Start with simple power-law, add complexity as needed |
| Performance (large N) | Barnes-Hut quadtree serves dual purpose (gravity + collisions) |
| Scope creep | Phases 1-4 are self-contained MVP |
| 2D limitations | Good enough for disk simulations, can extend to 3D later |

---

## Timeline Summary

| Week | Phase | Deliverable |
|------|-------|-------------|
| 0 | Core Infrastructure (Partial) | Body (2D), arena_bhtree quadtree ✅ |
| 1 | Core Infrastructure (Complete) | SystemState, orbital elements, tests |
| 2 | Gravity & Integration | Leapfrog, energy conservation |
| 3 | Collision Handling | Detection (direct, tree), merging, ejections |
| 4 | System I/O | Import/export with statistical crate |
| 5-6 | Gas Disk Effects | Drag, migration, damping |

**Recommended stopping points:**
- After Week 4: Functional pure-gravity N-body with I/O (tree-based collisions)
- After Week 6: Full implementation with gas effects

---

**Last Updated:** December 2025 (2D version)
