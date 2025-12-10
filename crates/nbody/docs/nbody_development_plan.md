# N-Body Engine Development Plan

## Overview

A gravitational N-body simulation engine that integrates with the existing Magrathea statistical and emergent crates. Designed for:

1. **Stability testing** — Evolve statistically-generated systems to verify dynamical plausibility
2. **Late-stage evolution** — Giant impact phase, resonant chain formation
3. **What-if experiments** — Perturb systems and observe outcomes
4. **Future gas effects** — Clean extension points for migration and drag

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
│           ├── body.rs           # Body, BodyId, OrbitalState
│           ├── state.rs          # SystemState, snapshot/restore
│           ├── vec3.rs           # Vector math (or use nalgebra)
│           ├── forces/
│           │   ├── mod.rs        # ForceModel trait, CompositeForce
│           │   ├── gravity.rs    # Direct N-body gravity
│           │   ├── tree.rs       # Barnes-Hut octree + neighbor queries
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

## Phase 1: Core Infrastructure (Week 1)

### Goals
- Basic data structures
- Vector math
- Body representation with orbital state

### 1.1 Vector Type

Either use `nalgebra` or roll a simple Vec3. For educational clarity and minimal dependencies, a simple implementation:

```rust
// src/vec3.rs

use std::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };
    
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    
    pub fn magnitude_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    
    pub fn magnitude(self) -> f64 {
        self.magnitude_squared().sqrt()
    }
    
    pub fn normalized(self) -> Self {
        let mag = self.magnitude();
        if mag > 0.0 {
            self / mag
        } else {
            Self::ZERO
        }
    }
    
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

// Implement Add, Sub, Mul<f64>, Div<f64>, Neg...
```

**Decision point:** If already using `nalgebra` elsewhere, use `Vector3<f64>` instead.

### 1.2 Body Representation

```rust
// src/body.rs

use crate::vec3::Vec3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyId(pub u32);

#[derive(Debug, Clone)]
pub struct Body {
    pub id: BodyId,
    pub mass: f64,              // kg
    pub radius: f64,            // m (physical radius for collisions)
    pub position: Vec3,         // m (heliocentric Cartesian)
    pub velocity: Vec3,         // m/s
}

impl Body {
    pub fn momentum(&self) -> Vec3 {
        self.velocity * self.mass
    }
    
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.magnitude_squared()
    }
    
    pub fn distance_to(&self, other: &Body) -> f64 {
        (self.position - other.position).magnitude()
    }
    
    pub fn orbital_radius(&self) -> f64 {
        self.position.magnitude()
    }
    
    pub fn orbital_velocity(&self) -> f64 {
        self.velocity.magnitude()
    }
    
    /// Angular momentum vector (r × v, not multiplied by mass)
    pub fn specific_angular_momentum(&self) -> Vec3 {
        self.position.cross(self.velocity)
    }
}
```

### 1.3 System State

```rust
// src/state.rs

use crate::body::{Body, BodyId};

#[derive(Debug, Clone)]
pub struct Star {
    pub mass: f64,      // kg
    pub radius: f64,    // m
}

#[derive(Debug, Clone)]
pub struct SystemState {
    pub time: f64,              // seconds since epoch
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
    
    pub fn add_body(&mut self, mass: f64, radius: f64, position: Vec3, velocity: Vec3) -> BodyId {
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

### 1.4 Orbital Elements Conversion

For diagnostics and I/O, convert between Cartesian and Keplerian:

```rust
// src/body.rs (continued)

use std::f64::consts::PI;

pub struct OrbitalElements {
    pub semi_major_axis: f64,   // m
    pub eccentricity: f64,
    pub inclination: f64,       // radians
    pub longitude_of_ascending_node: f64,  // radians
    pub argument_of_periapsis: f64,        // radians
    pub true_anomaly: f64,      // radians
}

impl Body {
    pub fn orbital_elements(&self, mu: f64) -> OrbitalElements {
        let r = self.position;
        let v = self.velocity;
        let r_mag = r.magnitude();
        let v_mag = v.magnitude();
        
        // Specific angular momentum
        let h = r.cross(v);
        let h_mag = h.magnitude();
        
        // Eccentricity vector
        let e_vec = v.cross(h) / mu - r.normalized();
        let e = e_vec.magnitude();
        
        // Semi-major axis (vis-viva)
        let energy = v_mag * v_mag / 2.0 - mu / r_mag;
        let a = -mu / (2.0 * energy);
        
        // Inclination
        let i = (h.z / h_mag).acos();
        
        // Node vector (z × h)
        let n = Vec3::new(-h.y, h.x, 0.0);
        let n_mag = n.magnitude();
        
        // Longitude of ascending node
        let omega = if n_mag > 1e-10 {
            let omega = (n.x / n_mag).acos();
            if n.y < 0.0 { 2.0 * PI - omega } else { omega }
        } else {
            0.0
        };
        
        // Argument of periapsis
        let w = if n_mag > 1e-10 && e > 1e-10 {
            let w = (n.dot(e_vec) / (n_mag * e)).acos();
            if e_vec.z < 0.0 { 2.0 * PI - w } else { w }
        } else {
            0.0
        };
        
        // True anomaly
        let nu = if e > 1e-10 {
            let nu = (e_vec.dot(r) / (e * r_mag)).clamp(-1.0, 1.0).acos();
            if r.dot(v) < 0.0 { 2.0 * PI - nu } else { nu }
        } else {
            0.0
        };
        
        OrbitalElements {
            semi_major_axis: a,
            eccentricity: e,
            inclination: i,
            longitude_of_ascending_node: omega,
            argument_of_periapsis: w,
            true_anomaly: nu,
        }
    }
}
```

### Phase 1 Deliverables

- [ ] `Vec3` with full operator overloads and geometric methods
- [ ] `Body` with momentum, energy, angular momentum
- [ ] `SystemState` with add/remove/lookup
- [ ] `OrbitalElements` conversion (Cartesian ↔ Keplerian)
- [ ] Unit tests for vector math
- [ ] Unit test: circular orbit elements should give e ≈ 0

### Phase 1 Test

```rust
#[test]
fn circular_orbit_elements() {
    let mu = G * M_SUN;
    let a = AU;  // 1 AU
    let v_circular = (mu / a).sqrt();
    
    let body = Body {
        id: BodyId(0),
        mass: M_EARTH,
        radius: R_EARTH,
        position: Vec3::new(a, 0.0, 0.0),
        velocity: Vec3::new(0.0, v_circular, 0.0),
    };
    
    let elements = body.orbital_elements(mu);
    
    assert!((elements.semi_major_axis - a).abs() / a < 1e-10);
    assert!(elements.eccentricity < 1e-10);
    assert!(elements.inclination.abs() < 1e-10);
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
use crate::vec3::Vec3;

/// A source of acceleration on bodies
pub trait ForceModel: Send + Sync {
    /// Compute acceleration on body at index `idx` given full system state
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vec3;
    
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
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vec3 {
        self.models.iter()
            .map(|f| f.acceleration(idx, state))
            .fold(Vec3::ZERO, |acc, a| acc + a)
    }
    
    fn potential_energy(&self, state: &SystemState) -> f64 {
        self.models.iter()
            .map(|f| f.potential_energy(state))
            .sum()
    }
}
```

### 2.2 Direct Gravity

```rust
// src/forces/gravity.rs

use crate::forces::ForceModel;
use crate::state::SystemState;
use crate::vec3::Vec3;

pub const G: f64 = 6.674e-11;  // m³ kg⁻¹ s⁻²

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
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vec3 {
        let body = &state.bodies[idx];
        let eps2 = self.softening * self.softening;
        
        // Acceleration from star
        let r_star = body.position;  // Star at origin
        let r2_star = r_star.magnitude_squared() + eps2;
        let r_star_mag = r2_star.sqrt();
        let a_star = -r_star * (G * state.star.mass / (r2_star * r_star_mag));
        
        // Acceleration from other bodies
        let a_bodies = state.bodies.iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .map(|(_, other)| {
                let dr = other.position - body.position;
                let r2 = dr.magnitude_squared() + eps2;
                let r = r2.sqrt();
                dr * (G * other.mass / (r2 * r))
            })
            .fold(Vec3::ZERO, |acc, a| acc + a);
        
        a_star + a_bodies
    }
    
    fn potential_energy(&self, state: &SystemState) -> f64 {
        let eps2 = self.softening * self.softening;
        
        // Star-body potential
        let star_potential: f64 = state.bodies.iter()
            .map(|b| {
                let r = (b.position.magnitude_squared() + eps2).sqrt();
                -G * state.star.mass * b.mass / r
            })
            .sum();
        
        // Body-body potential (each pair counted once)
        let body_potential: f64 = state.bodies.iter()
            .enumerate()
            .flat_map(|(i, a)| {
                state.bodies[i+1..].iter().map(move |b| {
                    let r = ((a.position - b.position).magnitude_squared() + eps2).sqrt();
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
use crate::vec3::Vec3;

pub struct Leapfrog<F: ForceModel> {
    force: F,
    dt: f64,
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
        let accelerations: Vec<Vec3> = (0..n)
            .map(|i| self.force.acceleration(i, state))
            .collect();
        
        // Kick: v += a * dt/2
        for (body, &acc) in state.bodies.iter_mut().zip(&accelerations) {
            body.velocity = body.velocity + acc * (dt / 2.0);
        }
        
        // Drift: r += v * dt
        for body in &mut state.bodies {
            body.position = body.position + body.velocity * dt;
        }
        
        // Recompute accelerations at new positions
        let accelerations: Vec<Vec3> = (0..n)
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
fn two_body_ellipse() {
    let mu = G * M_SUN;
    let a = AU;
    let e = 0.3;
    
    // Start at perihelion
    let r_peri = a * (1.0 - e);
    let v_peri = ((1.0 + e) * mu / (a * (1.0 - e))).sqrt();
    
    let star = Star { mass: M_SUN, radius: R_SUN };
    let mut state = SystemState::new(star);
    state.add_body(
        M_EARTH, R_EARTH,
        Vec3::new(r_peri, 0.0, 0.0),
        Vec3::new(0.0, v_peri, 0.0),
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
    let final_r = state.bodies[0].position.magnitude();
    assert!((final_r - r_peri).abs() / r_peri < 1e-4);
}
```

---

## Phase 3: Collision Handling (Week 3)

### Goals
- Detect close approaches
- Resolve collisions via merging
- Handle edge cases (ejections, star collisions)

### 3.1 Collision Detection

The naive O(N²) approach works for small N but becomes a bottleneck with thousands of bodies. We provide two spatial acceleration structures:

1. **Tree-based** (primary) — reuses the Barnes-Hut octree from gravity
2. **Radial bins** (lightweight) — exploits orbital geometry, simpler to implement

```rust
// src/collisions/detection.rs

use crate::body::{Body, BodyId};
use crate::state::SystemState;
use crate::vec3::Vec3;

pub const G: f64 = 6.674e-11;

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
// Radial Bins - Lightweight, Exploits Orbital Geometry
// ============================================================================

/// Logarithmically-spaced radial bins for broad-phase collision detection.
/// Effective for disk-like distributions where bodies are spread radially.
pub struct RadialBinDetector {
    /// Bin edges in meters (logarithmically spaced)
    bin_edges: Vec<f64>,
}

impl RadialBinDetector {
    /// Create with logarithmic bins from inner to outer edge
    pub fn new(inner_edge: f64, outer_edge: f64, bins_per_decade: usize) -> Self {
        let log_inner = inner_edge.log10();
        let log_outer = outer_edge.log10();
        let decades = log_outer - log_inner;
        let num_bins = (decades * bins_per_decade as f64).ceil() as usize;
        
        let bin_edges: Vec<f64> = (0..=num_bins)
            .map(|i| {
                let log_r = log_inner + (i as f64 / num_bins as f64) * (log_outer - log_inner);
                10f64.powf(log_r)
            })
            .collect();
        
        Self { bin_edges }
    }
    
    /// Default configuration suitable for planetary systems (0.1 AU to 1000 AU)
    pub fn default_planetary() -> Self {
        const AU: f64 = 1.496e11;
        Self::new(0.1 * AU, 1000.0 * AU, 10)  // 10 bins per decade
    }
    
    fn bin_index(&self, r: f64) -> Option<usize> {
        if r < self.bin_edges[0] || r >= *self.bin_edges.last().unwrap() {
            return None;
        }
        
        self.bin_edges
            .windows(2)
            .position(|w| r >= w[0] && r < w[1])
    }
    
    fn build_bins(&self, bodies: &[Body]) -> Vec<Vec<usize>> {
        let mut bins = vec![Vec::new(); self.bin_edges.len() - 1];
        
        for (idx, body) in bodies.iter().enumerate() {
            if let Some(bin) = self.bin_index(body.orbital_radius()) {
                bins[bin].push(idx);
            }
        }
        
        bins
    }
    
    /// Get bin indices that could contain colliders for a body at radius r with search_radius
    fn candidate_bins(&self, r: f64, search_radius: f64) -> std::ops::RangeInclusive<usize> {
        let r_min = (r - search_radius).max(self.bin_edges[0]);
        let r_max = r + search_radius;
        
        let bin_min = self.bin_index(r_min).unwrap_or(0);
        let bin_max = self.bin_index(r_max).unwrap_or(self.bin_edges.len() - 2);
        
        bin_min..=bin_max
    }
}

impl CollisionDetector for RadialBinDetector {
    fn detect(&self, state: &SystemState, criteria: &CollisionCriteria) -> Vec<CollisionEvent> {
        let bins = self.build_bins(&state.bodies);
        let mut events = Vec::new();
        let mut checked = std::collections::HashSet::new();
        
        for (idx, body) in state.bodies.iter().enumerate() {
            let r = body.orbital_radius();
            let r_hill = hill_radius(body.mass, r, state.star.mass);
            // Search radius: own Hill radius + maximum possible neighbor Hill radius
            // Conservative estimate: use 3x own Hill radius
            let search_radius = r_hill * 3.0;
            
            for bin_idx in self.candidate_bins(r, search_radius) {
                if let Some(bin) = bins.get(bin_idx) {
                    for &other_idx in bin {
                        if other_idx <= idx {
                            continue;  // Avoid self and duplicates
                        }
                        
                        let pair = (idx, other_idx);
                        if checked.contains(&pair) {
                            continue;
                        }
                        checked.insert(pair);
                        
                        if let Some(event) = check_pair(
                            &state.bodies[idx],
                            &state.bodies[other_idx],
                            state.star.mass,
                            criteria,
                        ) {
                            events.push(event);
                        }
                    }
                }
            }
        }
        
        events
    }
}

// ============================================================================
// Tree-Based - Reuses Barnes-Hut Octree
// ============================================================================

use crate::forces::tree::OctreeNode;

/// Collision detection using Barnes-Hut octree neighbor queries.
/// Most efficient when tree is already built for gravity calculation.
pub struct TreeDetector;

impl TreeDetector {
    /// Detect collisions using a pre-built octree
    pub fn detect_with_tree(
        &self,
        tree: &OctreeNode,
        state: &SystemState,
        criteria: &CollisionCriteria,
    ) -> Vec<CollisionEvent> {
        let mut events = Vec::new();
        let mut seen: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
        
        for body in &state.bodies {
            let r_hill = hill_radius(body.mass, body.orbital_radius(), state.star.mass);
            // Search for neighbors within ~2x Hill radius (to catch mutual Hill sphere)
            let search_radius = r_hill * 3.0;
            
            let mut neighbors = Vec::new();
            tree.query_radius(body.position, search_radius, &mut neighbors);
            
            for (neighbor_id, neighbor_pos) in neighbors {
                if neighbor_id == body.id {
                    continue;
                }
                
                // Canonical ordering to avoid duplicates
                let pair = if body.id.0 < neighbor_id.0 {
                    (body.id.0, neighbor_id.0)
                } else {
                    (neighbor_id.0, body.id.0)
                };
                
                if seen.contains(&pair) {
                    continue;
                }
                seen.insert(pair);
                
                // Get full body data for proper collision check
                if let Some(neighbor) = state.get_body(neighbor_id) {
                    if let Some(event) = check_pair(body, neighbor, state.star.mass, criteria) {
                        events.push(event);
                    }
                }
            }
        }
        
        events
    }
}

// ============================================================================
// Octree Neighbor Query Extension
// ============================================================================

// Add to src/forces/tree.rs:

impl OctreeNode {
    /// Find all bodies within `radius` of `position`
    pub fn query_radius(
        &self,
        position: Vec3,
        radius: f64,
        results: &mut Vec<(BodyId, Vec3)>,
    ) {
        match self {
            OctreeNode::Empty => {}
            
            OctreeNode::Leaf { body_id, position: leaf_pos, .. } => {
                let dist = (*leaf_pos - position).magnitude();
                if dist < radius {
                    results.push((*body_id, *leaf_pos));
                }
            }
            
            OctreeNode::Internal { children, bounds, .. } => {
                // Prune: skip if bounding box doesn't intersect search sphere
                if !bounds.intersects_sphere(position, radius) {
                    return;
                }
                
                for child in children.iter() {
                    child.query_radius(position, radius, results);
                }
            }
        }
    }
}

impl BoundingBox {
    /// Check if bounding box intersects a sphere
    pub fn intersects_sphere(&self, center: Vec3, radius: f64) -> bool {
        // Find closest point on box to sphere center
        let closest = Vec3::new(
            center.x.clamp(self.center.x - self.half_size, self.center.x + self.half_size),
            center.y.clamp(self.center.y - self.half_size, self.center.y + self.half_size),
            center.z.clamp(self.center.z - self.half_size, self.center.z + self.half_size),
        );
        
        (closest - center).magnitude_squared() < radius * radius
    }
}

// ============================================================================
// Detector Selection Helper
// ============================================================================

/// Choose appropriate detector based on body count
pub fn select_detector(body_count: usize, have_tree: bool) -> Box<dyn CollisionDetector> {
    match (body_count, have_tree) {
        (n, _) if n < 100 => Box::new(DirectDetector),
        (_, true) => Box::new(TreeDetector),  // Note: requires separate tree-based call
        (_, false) => Box::new(RadialBinDetector::default_planetary()),
    }
}

// ============================================================================
// Star Collision and Ejection Detection (unchanged)
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

### 3.2 Collision Resolution

```rust
// src/collisions/resolution.rs

use crate::body::{Body, BodyId};
use crate::state::SystemState;
use crate::vec3::Vec3;

/// Merge two bodies, conserving momentum
pub fn merge_bodies(a: &Body, b: &Body, new_id: BodyId) -> Body {
    let total_mass = a.mass + b.mass;
    
    // Center of mass position
    let position = (a.position * a.mass + b.position * b.mass) / total_mass;
    
    // Momentum-conserving velocity
    let velocity = (a.momentum() + b.momentum()) / total_mass;
    
    // Combined radius (volume-conserving, assuming same density)
    // V = (4/3)πr³, so r_new = (r_a³ + r_b³)^(1/3)
    let radius = (a.radius.powi(3) + b.radius.powi(3)).powf(1.0 / 3.0);
    
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
        let a = state.bodies.iter().find(|b| b.id == event.body_a).cloned();
        let b = state.bodies.iter().find(|b| b.id == event.body_b).cloned();
        
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

### 3.3 Integrated Evolution with Collisions

```rust
// src/integrators/mod.rs (addition)

use crate::collisions::{detection, resolution};
use crate::forces::ForceModel;
use crate::state::SystemState;

pub struct EvolutionConfig {
    pub collision_criteria: detection::CollisionCriteria,
    pub check_star_collision: bool,
    pub check_ejection: bool,
    pub ejection_radius: f64,  // typically 1000 AU
}

impl Default for EvolutionConfig {
    fn default() -> Self {
        Self {
            collision_criteria: detection::CollisionCriteria::default(),
            check_star_collision: true,
            check_ejection: true,
            ejection_radius: 1000.0 * 1.496e11,  // 1000 AU in meters
        }
    }
}

pub struct CollisionStats {
    pub mergers: u32,
    pub star_collisions: u32,
    pub ejections: u32,
}

pub fn evolve_with_collisions<F: ForceModel>(
    integrator: &super::leapfrog::Leapfrog<F>,
    state: &mut SystemState,
    duration: f64,
    config: &EvolutionConfig,
) -> CollisionStats {
    let mut stats = CollisionStats {
        mergers: 0,
        star_collisions: 0,
        ejections: 0,
    };
    
    let dt = integrator.dt;
    let steps = (duration / dt).ceil() as usize;
    
    for _ in 0..steps {
        integrator.step(state);
        
        // Check for collisions
        let collisions = detection::detect_collisions(state, &config.collision_criteria);
        stats.mergers += collisions.len() as u32;
        resolution::resolve_collisions(state, collisions);
        
        // Check for star collisions
        if config.check_star_collision {
            let star_hits = detection::detect_star_collisions(state);
            stats.star_collisions += star_hits.len() as u32;
            resolution::remove_star_collisions(state, star_hits);
        }
        
        // Check for ejections
        if config.check_ejection {
            let ejected = detection::detect_ejections(state, config.ejection_radius);
            stats.ejections += ejected.len() as u32;
            resolution::remove_ejections(state, ejected);
        }
    }
    
    stats
}
```

### Phase 3 Deliverables

- [ ] `CollisionDetector` trait with three implementations:
  - [ ] `DirectDetector` — O(N²) baseline for small N
  - [ ] `RadialBinDetector` — logarithmic bins, O(N × k) 
  - [ ] `TreeDetector` — reuses Barnes-Hut octree
- [ ] `BoundingBox::intersects_sphere()` for tree queries
- [ ] `OctreeNode::query_radius()` for neighbor search
- [ ] `CollisionCriteria` configuration
- [ ] `merge_bodies()` with momentum conservation
- [ ] `detect_star_collisions()` and `detect_ejections()`
- [ ] `evolve_with_collisions()` integrated loop
- [ ] Test: Two bodies on collision course merge correctly
- [ ] Test: Body on escape trajectory is removed
- [ ] Benchmark: Measure crossover point where radial bins beat direct

### Phase 3 Test: Collision and Merger

```rust
#[test]
fn collision_conserves_momentum() {
    let star = Star { mass: M_SUN, radius: R_SUN };
    let mut state = SystemState::new(star);
    
    // Two bodies approaching each other
    state.add_body(
        M_EARTH, R_EARTH,
        Vec3::new(AU, 0.0, 0.0),
        Vec3::new(0.0, 29780.0, 0.0),  // ~circular at 1 AU
    );
    state.add_body(
        M_EARTH, R_EARTH,
        Vec3::new(AU + 2.0 * R_EARTH, 0.0, 0.0),  // Very close
        Vec3::new(0.0, 29780.0, 0.0),
    );
    
    let initial_momentum = state.bodies.iter()
        .map(|b| b.momentum())
        .fold(Vec3::ZERO, |a, b| a + b);
    
    // Detect and resolve
    let events = detection::detect_collisions(&state, &detection::CollisionCriteria::default());
    assert_eq!(events.len(), 1);
    
    resolution::resolve_collisions(&mut state, events);
    assert_eq!(state.bodies.len(), 1);
    
    let final_momentum = state.bodies[0].momentum();
    
    // Momentum conserved
    assert!((final_momentum - initial_momentum).magnitude() < 1e-6);
    
    // Mass conserved
    assert!((state.bodies[0].mass - 2.0 * M_EARTH).abs() < 1e-6);
}
```

---

## Phase 4: System I/O (Week 4)

### Goals
- Import from existing `System` JSON format
- Export evolved state to `System` JSON
- Round-trip validation

### 4.1 Import from System

```rust
// src/io/import.rs

use crate::body::Body;
use crate::state::{Star, SystemState};
use crate::vec3::Vec3;
use std::f64::consts::PI;

// Assuming System is available from the planetary crate
// These would be the actual imports in practice:
// use planetary::System;
// use stellar::Star as StellarStar;

pub const G: f64 = 6.674e-11;
pub const AU: f64 = 1.496e11;

/// Convert orbital elements to Cartesian state
/// Assumes orbits in the x-y plane for simplicity (can extend to full 3D)
pub fn elements_to_cartesian(
    a: f64,           // semi-major axis (m)
    e: f64,           // eccentricity
    i: f64,           // inclination (rad)
    omega: f64,       // longitude of ascending node (rad)
    w: f64,           // argument of periapsis (rad)
    nu: f64,          // true anomaly (rad)
    mu: f64,          // G * (M_star + m_body)
) -> (Vec3, Vec3) {
    // Distance from focus
    let r = a * (1.0 - e * e) / (1.0 + e * nu.cos());
    
    // Position and velocity in orbital plane
    let x_orb = r * nu.cos();
    let y_orb = r * nu.sin();
    
    let h = (mu * a * (1.0 - e * e)).sqrt();  // specific angular momentum
    let vx_orb = -mu / h * nu.sin();
    let vy_orb = mu / h * (e + nu.cos());
    
    // Rotation matrices (combined)
    let cos_o = omega.cos();
    let sin_o = omega.sin();
    let cos_w = w.cos();
    let sin_w = w.sin();
    let cos_i = i.cos();
    let sin_i = i.sin();
    
    // Transform to inertial frame
    let x = (cos_o * cos_w - sin_o * sin_w * cos_i) * x_orb 
          + (-cos_o * sin_w - sin_o * cos_w * cos_i) * y_orb;
    let y = (sin_o * cos_w + cos_o * sin_w * cos_i) * x_orb 
          + (-sin_o * sin_w + cos_o * cos_w * cos_i) * y_orb;
    let z = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb;
    
    let vx = (cos_o * cos_w - sin_o * sin_w * cos_i) * vx_orb 
           + (-cos_o * sin_w - sin_o * cos_w * cos_i) * vy_orb;
    let vy = (sin_o * cos_w + cos_o * sin_w * cos_i) * vx_orb 
           + (-sin_o * sin_w + cos_o * cos_w * cos_i) * vy_orb;
    let vz = (sin_w * sin_i) * vx_orb + (cos_w * sin_i) * vy_orb;
    
    (Vec3::new(x, y, z), Vec3::new(vx, vy, vz))
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
    pub inclination: f64,     // degrees
}

pub const M_SUN: f64 = 1.989e30;
pub const R_SUN: f64 = 6.957e8;
pub const M_EARTH: f64 = 5.972e24;
pub const R_EARTH: f64 = 6.371e6;

pub fn import_system(json: &SystemJson, config: &ImportConfig) -> SystemState {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    
    let star = Star {
        mass: json.star_mass * M_SUN,
        radius: json.star_radius * R_SUN,
    };
    
    let mut state = SystemState::new(star);
    let mut rng = ChaCha8Rng::seed_from_u64(config.seed);
    
    for (i, planet) in json.planets.iter().enumerate() {
        let a = planet.semi_major_axis * AU;
        let e = planet.eccentricity;
        let inc = planet.inclination.to_radians();
        let mass = planet.mass * M_EARTH;
        let radius = planet.radius * R_EARTH;
        
        // Determine true anomaly
        let nu = match config.anomaly_mode {
            AnomalyMode::Periapsis => 0.0,
            AnomalyMode::Apoapsis => PI,
            AnomalyMode::Random => rng.gen::<f64>() * 2.0 * PI,
            AnomalyMode::Uniform => (i as f64 / json.planets.len() as f64) * 2.0 * PI,
        };
        
        // Random orbital angles (could also be imported if available)
        let omega = rng.gen::<f64>() * 2.0 * PI;
        let w = rng.gen::<f64>() * 2.0 * PI;
        
        let mu = G * (star.mass + mass);
        let (position, velocity) = elements_to_cartesian(a, e, inc, omega, w, nu, mu);
        
        state.add_body(mass, radius, position, velocity);
    }
    
    state
}
```

### 4.2 Export to System

```rust
// src/io/export.rs

use crate::state::SystemState;
use crate::body::OrbitalElements;

pub const G: f64 = 6.674e-11;
pub const M_SUN: f64 = 1.989e30;
pub const R_SUN: f64 = 6.957e8;
pub const M_EARTH: f64 = 5.972e24;
pub const R_EARTH: f64 = 6.371e6;
pub const AU: f64 = 1.496e11;

use super::import::{SystemJson, PlanetJson};

pub fn export_system(state: &SystemState) -> SystemJson {
    let mu = G * state.star.mass;
    
    let planets = state.bodies.iter()
        .map(|body| {
            let elements = body.orbital_elements(mu);
            
            PlanetJson {
                mass: body.mass / M_EARTH,
                radius: body.radius / R_EARTH,
                semi_major_axis: elements.semi_major_axis / AU,
                eccentricity: elements.eccentricity,
                inclination: elements.inclination.to_degrees(),
            }
        })
        .collect();
    
    SystemJson {
        star_mass: state.star.mass / M_SUN,
        star_radius: state.star.radius / R_SUN,
        planets,
    }
}
```

### Phase 4 Deliverables

- [ ] `elements_to_cartesian()` conversion
- [ ] `import_system()` with configurable anomaly modes
- [ ] `export_system()` to JSON
- [ ] Round-trip test: import → export preserves orbital elements
- [ ] Integration test: import statistical output, evolve, export

---

## Phase 5: Gas Disk Effects (Week 5-6)

### Goals
- Pluggable gas disk models
- Aerodynamic drag for small bodies
- Type I migration for planets
- Eccentricity/inclination damping

### 5.1 Gas Disk Trait

```rust
// src/forces/gas/disk.rs

/// Properties of a gas disk at a given location and time
pub trait GasDiskModel: Send + Sync {
    /// Surface density (kg/m²)
    fn surface_density(&self, r: f64, t: f64) -> f64;
    
    /// Midplane temperature (K)
    fn temperature(&self, r: f64, t: f64) -> f64;
    
    /// Scale height (m)
    fn scale_height(&self, r: f64, t: f64) -> f64;
    
    /// Aspect ratio h/r
    fn aspect_ratio(&self, r: f64, t: f64) -> f64 {
        self.scale_height(r, t) / r
    }
    
    /// Midplane gas density (kg/m³)
    fn midplane_density(&self, r: f64, t: f64) -> f64 {
        // ρ = Σ / (√(2π) × H)
        self.surface_density(r, t) / (2.506628 * self.scale_height(r, t))
    }
    
    /// Gas velocity (slightly sub-Keplerian due to pressure support)
    fn gas_velocity(&self, r: f64, t: f64, stellar_mass: f64) -> f64;
    
    /// Is the disk still present? (for early termination)
    fn is_dispersed(&self, t: f64) -> bool;
}
```

### 5.2 Power-Law Disk Implementation

```rust
// src/forces/gas/disk.rs (continued)

pub const K_BOLTZMANN: f64 = 1.381e-23;
pub const M_PROTON: f64 = 1.673e-27;
pub const MU_GAS: f64 = 2.34;  // Mean molecular weight
pub const AU: f64 = 1.496e11;
pub const G: f64 = 6.674e-11;

/// Simple power-law disk with exponential decay
pub struct PowerLawDisk {
    /// Surface density at 1 AU (kg/m²)
    pub sigma_0: f64,
    /// Surface density power law exponent (typically -1 to -1.5)
    pub sigma_exp: f64,
    /// Temperature at 1 AU (K)
    pub temp_0: f64,
    /// Temperature power law exponent (typically -0.5)
    pub temp_exp: f64,
    /// Disk dispersal timescale (s)
    pub dispersal_time: f64,
    /// Inner edge (m)
    pub inner_edge: f64,
    /// Outer edge (m)
    pub outer_edge: f64,
    /// Stellar mass (kg) - needed for velocity calculation
    pub stellar_mass: f64,
}

impl PowerLawDisk {
    /// MMSN-like disk around a solar-mass star
    pub fn mmsn(stellar_mass: f64) -> Self {
        Self {
            sigma_0: 1700.0,           // kg/m² at 1 AU
            sigma_exp: -1.5,
            temp_0: 280.0,             // K at 1 AU
            temp_exp: -0.5,
            dispersal_time: 3.0e6 * 3.156e7,  // 3 Myr in seconds
            inner_edge: 0.1 * AU,
            outer_edge: 100.0 * AU,
            stellar_mass,
        }
    }
    
    fn sound_speed(&self, r: f64, t: f64) -> f64 {
        let temp = self.temperature(r, t);
        (K_BOLTZMANN * temp / (MU_GAS * M_PROTON)).sqrt()
    }
    
    fn pressure_gradient_eta(&self, r: f64, t: f64) -> f64 {
        // η = -(h/r)² × (d ln P / d ln r) / 2
        // For power-law disk: d ln P / d ln r ≈ sigma_exp + temp_exp - 1.5
        let h_r = self.aspect_ratio(r, t);
        let p_exp = self.sigma_exp + self.temp_exp - 1.5;
        -h_r * h_r * p_exp / 2.0
    }
}

impl GasDiskModel for PowerLawDisk {
    fn surface_density(&self, r: f64, t: f64) -> f64 {
        if r < self.inner_edge || r > self.outer_edge {
            return 0.0;
        }
        let decay = (-t / self.dispersal_time).exp();
        self.sigma_0 * (r / AU).powf(self.sigma_exp) * decay
    }
    
    fn temperature(&self, r: f64, _t: f64) -> f64 {
        // Temperature doesn't decay (set by stellar irradiation)
        self.temp_0 * (r / AU).powf(self.temp_exp)
    }
    
    fn scale_height(&self, r: f64, t: f64) -> f64 {
        let cs = self.sound_speed(r, t);
        let omega = (G * self.stellar_mass / r.powi(3)).sqrt();
        cs / omega
    }
    
    fn gas_velocity(&self, r: f64, t: f64, stellar_mass: f64) -> f64 {
        let v_kep = (G * stellar_mass / r).sqrt();
        let eta = self.pressure_gradient_eta(r, t);
        v_kep * (1.0 - eta)
    }
    
    fn is_dispersed(&self, t: f64) -> bool {
        // Consider dispersed after 5 e-folding times
        t > 5.0 * self.dispersal_time
    }
}
```

### 5.3 Drag Force

```rust
// src/forces/gas/drag.rs

use crate::forces::ForceModel;
use crate::state::SystemState;
use crate::vec3::Vec3;
use super::disk::GasDiskModel;

pub struct AerodynamicDrag<D: GasDiskModel> {
    pub disk: D,
    /// Drag coefficient (typically ~1)
    pub c_d: f64,
    /// Maximum mass affected by drag (kg)
    pub max_mass: f64,
}

impl<D: GasDiskModel> AerodynamicDrag<D> {
    pub fn new(disk: D) -> Self {
        Self {
            disk,
            c_d: 1.0,
            max_mass: 1e21,  // ~100 km body
        }
    }
}

impl<D: GasDiskModel> ForceModel for AerodynamicDrag<D> {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vec3 {
        let body = &state.bodies[idx];
        
        // Skip massive bodies
        if body.mass > self.max_mass {
            return Vec3::ZERO;
        }
        
        let r = body.orbital_radius();
        let t = state.time;
        
        let rho_gas = self.disk.midplane_density(r, t);
        if rho_gas <= 0.0 {
            return Vec3::ZERO;
        }
        
        // Gas velocity (circular, in x-y plane)
        let v_gas_mag = self.disk.gas_velocity(r, t, state.star.mass);
        let r_hat = body.position.normalized();
        let z_hat = Vec3::new(0.0, 0.0, 1.0);
        let phi_hat = z_hat.cross(r_hat).normalized();
        let v_gas = phi_hat * v_gas_mag;
        
        // Relative velocity
        let delta_v = body.velocity - v_gas;
        let delta_v_mag = delta_v.magnitude();
        
        if delta_v_mag < 1e-10 {
            return Vec3::ZERO;
        }
        
        // Cross-sectional area
        let area = std::f64::consts::PI * body.radius * body.radius;
        
        // Drag acceleration: a = -½ × C_D × ρ × A × |Δv| × Δv / m
        let drag_mag = 0.5 * self.c_d * rho_gas * area * delta_v_mag / body.mass;
        
        -delta_v.normalized() * drag_mag * delta_v_mag
    }
}
```

### 5.4 Type I Migration

```rust
// src/forces/gas/migration.rs

use crate::forces::ForceModel;
use crate::state::SystemState;
use crate::vec3::Vec3;
use super::disk::GasDiskModel;

pub const G: f64 = 6.674e-11;

pub struct TypeIMigration<D: GasDiskModel> {
    pub disk: D,
    /// Minimum mass for migration (kg)
    pub min_mass: f64,
    /// Migration efficiency factor (0-1, for tuning)
    pub efficiency: f64,
}

impl<D: GasDiskModel> TypeIMigration<D> {
    pub fn new(disk: D) -> Self {
        Self {
            disk,
            min_mass: 1e23,  // ~Moon mass
            efficiency: 1.0,
        }
    }
    
    /// Migration timescale (Tanaka et al. 2002)
    fn migration_timescale(&self, body_mass: f64, r: f64, t: f64, stellar_mass: f64) -> f64 {
        let sigma = self.disk.surface_density(r, t);
        if sigma <= 0.0 {
            return f64::INFINITY;
        }
        
        let h_r = self.disk.aspect_ratio(r, t);
        let omega = (G * stellar_mass / r.powi(3)).sqrt();
        
        // τ = (M*/m) × (M*/(Σr²)) × (h/r)² × Ω⁻¹
        let m_ratio = stellar_mass / body_mass;
        let sigma_factor = stellar_mass / (sigma * r * r);
        
        m_ratio * sigma_factor * h_r * h_r / omega / self.efficiency
    }
    
    /// Eccentricity damping timescale (typically τ_e ~ τ_mig / 10)
    fn eccentricity_damping_timescale(&self, body_mass: f64, r: f64, t: f64, stellar_mass: f64) -> f64 {
        self.migration_timescale(body_mass, r, t, stellar_mass) / 10.0
    }
    
    /// Inclination damping timescale (typically τ_i ~ τ_e)
    fn inclination_damping_timescale(&self, body_mass: f64, r: f64, t: f64, stellar_mass: f64) -> f64 {
        self.eccentricity_damping_timescale(body_mass, r, t, stellar_mass)
    }
}

impl<D: GasDiskModel> ForceModel for TypeIMigration<D> {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vec3 {
        let body = &state.bodies[idx];
        
        // Skip low-mass bodies (handled by drag instead)
        if body.mass < self.min_mass {
            return Vec3::ZERO;
        }
        
        let r = body.orbital_radius();
        let t = state.time;
        
        // Check if disk is present
        if self.disk.surface_density(r, t) <= 0.0 {
            return Vec3::ZERO;
        }
        
        // Decompose velocity into components
        let r_hat = body.position.normalized();
        let z_hat = Vec3::new(0.0, 0.0, 1.0);
        let phi_hat = z_hat.cross(r_hat).normalized();
        
        let v_r = body.velocity.dot(r_hat);      // Radial
        let v_phi = body.velocity.dot(phi_hat);  // Tangential
        let v_z = body.velocity.dot(z_hat);      // Vertical
        
        // Timescales
        let tau_mig = self.migration_timescale(body.mass, r, t, state.star.mass);
        let tau_e = self.eccentricity_damping_timescale(body.mass, r, t, state.star.mass);
        let tau_i = self.inclination_damping_timescale(body.mass, r, t, state.star.mass);
        
        // Migration: damp tangential velocity slightly (causes inward drift)
        // Actually: apply a gentle inward acceleration
        let a_mig = if tau_mig.is_finite() {
            -r_hat * (r / tau_mig / tau_mig).sqrt()  // Approximate inward acceleration
        } else {
            Vec3::ZERO
        };
        
        // Eccentricity damping: damp radial velocity
        let a_e = if tau_e.is_finite() {
            -r_hat * (v_r / tau_e)
        } else {
            Vec3::ZERO
        };
        
        // Inclination damping: damp vertical velocity
        let a_i = if tau_i.is_finite() {
            -z_hat * (v_z / tau_i)
        } else {
            Vec3::ZERO
        };
        
        a_mig + a_e + a_i
    }
}
```

### Phase 5 Deliverables

- [ ] `GasDiskModel` trait
- [ ] `PowerLawDisk` implementation with MMSN defaults
- [ ] `AerodynamicDrag` force model for small bodies
- [ ] `TypeIMigration` force model with e/i damping
- [ ] Test: Planetesimal spirals inward due to drag
- [ ] Test: Earth-mass planet migrates inward
- [ ] Test: Disk dispersal stops migration

---

## Phase 6: Barnes-Hut Tree (Week 7-8)

### Goals
- Octree construction
- Tree-walk force calculation  
- Integrate with collision detection (already using tree queries)

The tree now serves dual purposes: O(N log N) gravity calculation AND efficient neighbor queries for collision detection. This phase focuses on the gravity side.

### 6.1 Octree Node

```rust
// src/forces/tree.rs

use crate::vec3::Vec3;
use crate::body::BodyId;

/// Bounding box in 3D
#[derive(Debug, Clone)]
pub struct BoundingBox {
    pub center: Vec3,
    pub half_size: f64,
}

impl BoundingBox {
    pub fn contains(&self, point: Vec3) -> bool {
        (point.x - self.center.x).abs() <= self.half_size
            && (point.y - self.center.y).abs() <= self.half_size
            && (point.z - self.center.z).abs() <= self.half_size
    }
    
    pub fn octant(&self, point: Vec3) -> usize {
        let mut index = 0;
        if point.x > self.center.x { index |= 1; }
        if point.y > self.center.y { index |= 2; }
        if point.z > self.center.z { index |= 4; }
        index
    }
    
    pub fn child_box(&self, octant: usize) -> BoundingBox {
        let offset = self.half_size / 2.0;
        let dx = if octant & 1 != 0 { offset } else { -offset };
        let dy = if octant & 2 != 0 { offset } else { -offset };
        let dz = if octant & 4 != 0 { offset } else { -offset };
        
        BoundingBox {
            center: Vec3::new(
                self.center.x + dx,
                self.center.y + dy,
                self.center.z + dz,
            ),
            half_size: offset,
        }
    }
    
    /// Check if bounding box intersects a sphere (for collision queries)
    pub fn intersects_sphere(&self, center: Vec3, radius: f64) -> bool {
        let closest = Vec3::new(
            center.x.clamp(self.center.x - self.half_size, self.center.x + self.half_size),
            center.y.clamp(self.center.y - self.half_size, self.center.y + self.half_size),
            center.z.clamp(self.center.z - self.half_size, self.center.z + self.half_size),
        );
        
        (closest - center).magnitude_squared() < radius * radius
    }
}

/// Octree node for Barnes-Hut
pub enum OctreeNode {
    Empty,
    Leaf {
        body_id: BodyId,
        position: Vec3,
        mass: f64,
    },
    Internal {
        children: Box<[OctreeNode; 8]>,
        center_of_mass: Vec3,
        total_mass: f64,
        bounds: BoundingBox,
    },
}

impl OctreeNode {
    pub fn new() -> Self {
        OctreeNode::Empty
    }
    
    // Insert, build, and walk methods...
}
```

### 6.2 Tree Construction

```rust
impl OctreeNode {
    pub fn build(bodies: &[(BodyId, Vec3, f64)], bounds: BoundingBox) -> Self {
        match bodies.len() {
            0 => OctreeNode::Empty,
            1 => OctreeNode::Leaf {
                body_id: bodies[0].0,
                position: bodies[0].1,
                mass: bodies[0].2,
            },
            _ => {
                // Partition bodies into octants
                let mut octant_bodies: [Vec<(BodyId, Vec3, f64)>; 8] = Default::default();
                
                for &(id, pos, mass) in bodies {
                    let octant = bounds.octant(pos);
                    octant_bodies[octant].push((id, pos, mass));
                }
                
                // Build children recursively
                let children: [OctreeNode; 8] = std::array::from_fn(|i| {
                    let child_bounds = bounds.child_box(i);
                    OctreeNode::build(&octant_bodies[i], child_bounds)
                });
                
                // Compute center of mass
                let total_mass: f64 = bodies.iter().map(|(_, _, m)| m).sum();
                let center_of_mass = bodies.iter()
                    .map(|(_, pos, m)| *pos * *m)
                    .fold(Vec3::ZERO, |a, b| a + b) / total_mass;
                
                OctreeNode::Internal {
                    children: Box::new(children),
                    center_of_mass,
                    total_mass,
                    bounds,
                }
            }
        }
    }
}
```

### 6.3 Tree Walk for Gravity

```rust
const G: f64 = 6.674e-11;

impl OctreeNode {
    /// Compute acceleration on a body using Barnes-Hut approximation
    pub fn acceleration(&self, position: Vec3, theta: f64, softening: f64) -> Vec3 {
        match self {
            OctreeNode::Empty => Vec3::ZERO,
            
            OctreeNode::Leaf { position: leaf_pos, mass, .. } => {
                let dr = *leaf_pos - position;
                let r2 = dr.magnitude_squared() + softening * softening;
                let r = r2.sqrt();
                
                if r < 1e-10 {
                    return Vec3::ZERO;  // Self-interaction
                }
                
                dr * (G * mass / (r2 * r))
            }
            
            OctreeNode::Internal { children, center_of_mass, total_mass, bounds } => {
                let dr = *center_of_mass - position;
                let r = dr.magnitude();
                
                // Opening criterion: size / distance < theta
                if bounds.half_size * 2.0 / r < theta {
                    // Use monopole approximation
                    let r2 = r * r + softening * softening;
                    dr * (G * total_mass / (r2 * r2.sqrt()))
                } else {
                    // Recurse into children
                    children.iter()
                        .map(|child| child.acceleration(position, theta, softening))
                        .fold(Vec3::ZERO, |a, b| a + b)
                }
            }
        }
    }
    
    /// Find all bodies within `radius` of `position` (for collision detection)
    pub fn query_radius(
        &self,
        position: Vec3,
        radius: f64,
        results: &mut Vec<(BodyId, Vec3)>,
    ) {
        match self {
            OctreeNode::Empty => {}
            
            OctreeNode::Leaf { body_id, position: leaf_pos, .. } => {
                let dist = (*leaf_pos - position).magnitude();
                if dist < radius {
                    results.push((*body_id, *leaf_pos));
                }
            }
            
            OctreeNode::Internal { children, bounds, .. } => {
                if !bounds.intersects_sphere(position, radius) {
                    return;
                }
                
                for child in children.iter() {
                    child.query_radius(position, radius, results);
                }
            }
        }
    }
}
```

### 6.4 Barnes-Hut Force Model

```rust
pub struct BarnesHutGravity {
    /// Opening angle (typically 0.5 - 1.0)
    pub theta: f64,
    /// Softening length
    pub softening: f64,
}

impl BarnesHutGravity {
    pub fn new(theta: f64) -> Self {
        Self { theta, softening: 0.0 }
    }
    
    pub fn with_softening(theta: f64, softening: f64) -> Self {
        Self { theta, softening }
    }
    
    /// Build tree and compute all accelerations
    pub fn accelerations(&self, state: &SystemState) -> (OctreeNode, Vec<Vec3>) {
        // Determine bounding box
        let max_r = state.bodies.iter()
            .map(|b| b.position.magnitude())
            .fold(0.0, f64::max);
        
        let bounds = BoundingBox {
            center: Vec3::ZERO,
            half_size: max_r * 1.1,
        };
        
        // Build tree
        let body_data: Vec<_> = state.bodies.iter()
            .map(|b| (b.id, b.position, b.mass))
            .collect();
        let tree = OctreeNode::build(&body_data, bounds);
        
        // Compute accelerations (tree + star)
        let mu_star = G * state.star.mass;
        let accelerations = state.bodies.iter()
            .map(|body| {
                // Star contribution
                let r = body.position;
                let r2 = r.magnitude_squared() + self.softening * self.softening;
                let r_mag = r2.sqrt();
                let a_star = -r * (mu_star / (r2 * r_mag));
                
                // Body contributions via tree
                let a_bodies = tree.acceleration(body.position, self.theta, self.softening);
                
                a_star + a_bodies
            })
            .collect();
        
        (tree, accelerations)
    }
}
```

### Phase 6 Deliverables

- [ ] `BoundingBox` with octant subdivision
- [ ] `OctreeNode` enum with Empty/Leaf/Internal
- [ ] Tree construction from body list
- [ ] Tree-walk acceleration with opening angle θ
- [ ] `BarnesHutGravity` implementing `ForceModel`
- [ ] Benchmark: Find crossover N where tree beats direct

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
- REBOUND (established N-body code)
- Known Solar System ephemerides
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
- [ ] Can simulate 10,000 planetesimals (with tree)
- [ ] Reproduce known migration timescales within factor of 2

---

## Dependencies

### Required
- `rand` + `rand_chacha` for seeded RNG (likely already have)

### Optional
- `nalgebra` for optimized vector math (if not rolling own)
- `rayon` for parallel force calculation
- `serde` for JSON I/O (likely already have)

### Development
- `criterion` for benchmarking
- `proptest` for property-based testing

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Integration instability | Use symplectic integrator, verify energy conservation |
| Close encounter singularities | Collision detection + merging before singularity |
| Gas physics too complex | Start with simple power-law, add complexity as needed |
| Performance (large N) | Barnes-Hut tree serves dual purpose (gravity + collisions) |
| Scope creep | Phases 1-4 are self-contained MVP |

---

## Timeline Summary

| Week | Phase | Deliverable |
|------|-------|-------------|
| 1 | Core Infrastructure | Bodies, vectors, orbital elements |
| 2 | Gravity & Integration | Leapfrog, energy conservation |
| 3 | Collision Handling | Detection (direct, radial bins), merging, ejections |
| 4 | System I/O | Import/export with statistical crate |
| 5-6 | Gas Disk Effects | Drag, migration, damping |
| 7-8 | Barnes-Hut Tree | O(N log N) gravity + tree-based collision detection |

**Recommended stopping points:**
- After Week 4: Functional pure-gravity N-body with I/O (radial bins for collisions)
- After Week 6: Full implementation with gas effects
- After Week 8: High-performance version for large N (tree for both gravity and collisions)
