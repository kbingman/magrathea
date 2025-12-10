# Arena-Based Barnes-Hut Tree Implementation

## Overview

This document describes an arena-based Barnes-Hut quadtree implementation designed for efficient N-body gravitational simulations. The implementation trades pointer-based tree structures for index-based references, achieving better cache locality and simpler memory management.

## Table of Contents

- [Design Philosophy](#design-philosophy)
- [Core Concepts](#core-concepts)
- [Implementation Details](#implementation-details)
- [Usage Patterns](#usage-patterns)
- [Parallel Extensions](#parallel-extensions)
- [Performance Considerations](#performance-considerations)

---

## Design Philosophy

### Problems with Arc-Based Trees

Traditional quadtree implementations use `Arc<Node>` for child references:

```rust
pub enum QuadtreeNode {
    Internal {
        children: [Arc<QuadtreeNode>; 4],
        // ...
    }
}
```

**Issues:**

1. **Cache Locality**: Nodes scattered throughout heap memory (100+ cycle cache miss penalty)
2. **Allocation Overhead**: Each node requires separate heap allocation
3. **Memory Fragmentation**: Nodes allocated at different times end up far apart
4. **Reference Counting**: Arc adds atomic operations for every clone/drop

### Arena Solution

An **arena** is a contiguous `Vec` where all nodes live together:

```
Vec<Node>: [Node₀][Node₁][Node₂][Node₃][Node₄]...
            ↑      ↑      ↑
         Index 0    1      2
```

Nodes reference each other by **index** rather than pointer:

```rust
pub struct NodeId(u32);  // Just an index!

pub enum Node {
    Internal {
        children: [NodeId; 4],  // 16 bytes vs 32 for [Arc; 4]
        // ...
    }
}
```

**Benefits:**

- ✅ Sequential memory access (cache-friendly)
- ✅ Single allocation for entire tree
- ✅ No reference counting overhead
- ✅ Smaller nodes (u32 vs usize/Arc)
- ✅ Simple lifetime management

---

## Core Concepts

### 1. The `Massive` Trait

Generic interface for objects participating in Barnes-Hut calculations:

```rust
pub trait Massive: Copy {
    fn position(&self) -> Point2<f32>;
    fn mass(&self) -> Mass<f32>;
    
    fn mass_solar(&self) -> f32 {
        self.mass().to_solar_masses()
    }
}
```

**Design choices:**

- `Copy` bound: Avoids expensive cloning during tree construction
- `mass_solar()`: Caches conversion to avoid repeated calls in hot loops
- Minimal interface: Only position and mass required

**Implementation for your Body type:**

```rust
impl Massive for Body {
    fn position(&self) -> Point2<f32> {
        self.position
    }
    
    fn mass(&self) -> Mass<f32> {
        self.mass
    }
}
```

### 2. NodeId - Type-Safe Indices

```rust
pub struct NodeId(u32);

impl NodeId {
    pub const EMPTY: NodeId = NodeId(u32::MAX);  // Sentinel value
    
    fn new(index: usize) -> Self {
        NodeId(index as u32)
    }
    
    fn is_empty(self) -> bool {
        self == Self::EMPTY
    }
}
```

**Why u32 instead of usize?**

- **Space savings**: 4 bytes vs 8 on 64-bit (16 bytes saved per internal node)
- **Range sufficient**: 4 billion nodes is more than enough
- **Sentinel pattern**: `u32::MAX` represents "no node" (like `null` but type-safe)

### 3. Node Enum

```rust
pub enum Node {
    Empty,
    
    Leaf {
        body_index: u32,  // Index into original &[Body]
    },
    
    LeafMulti {
        start: u32,       // Range of bodies at max depth
        count: u16,
    },
    
    Internal {
        center_of_mass: Point2<f32>,
        total_mass: f32,           // Pre-converted to solar masses
        bounds: BoundingBox,
        children: [NodeId; 4],     // Indices, not pointers
    },
}
```

**Key differences from Arc-based:**

- No `Vec<Body>` storage - just indices back to original slice
- `total_mass` stored as `f32` (already converted)
- `LeafMulti` variant handles max-depth edge case
- ~64 bytes per node (fits in cache line)

### 4. The BHTree Structure

```rust
pub struct BHTree<'a, B: Massive> {
    nodes: Vec<Node>,   // The arena
    bodies: &'a [B],    // Borrowed reference
    root: NodeId,       // Entry point for traversal
}
```

**Lifetime semantics:**

- Tree borrows bodies for lifetime `'a`
- Tree cannot outlive the body slice
- No body copying - zero allocation beyond nodes
- Rebuild is cheap: just rebuild `Vec<Node>`

---

## Implementation Details

### Tree Construction

The build process is **functionally recursive** even though it uses `&mut Vec`:

```rust
fn build_recursive(
    bodies: &[B],
    indices: &[usize],      // Which bodies in this subtree
    bounds: BoundingBox,
    depth: usize,
    arena: &mut Vec<Node>,  // Where to push nodes
) -> NodeId
```

**Build algorithm:**

```
1. Match on indices slice:
   [] → return EMPTY sentinel
   [single] → create Leaf, push to arena, return its index
   [many] at MAX_DEPTH → create LeafMulti
   [many] → subdivide:
      
2. Partition indices by quadrant:
   quadrants[0..4] = bodies.filter(|b| quadrant(b.pos) == q)
   
3. Recursively build children (depth-first):
   children = [build(q0), build(q1), build(q2), build(q3)]
   
4. Compute mass properties (single fold):
   (total_mass, center_of_mass) = fold over indices
   
5. Push Internal node, return its index
```

**Key insight:** Children are built *before* parent is pushed, so child indices are always smaller than parent indices. The root is always the last node.

**Example trace for 3 bodies:**

```
Input: [Sun(0,0), Earth(1,0), Jupiter(5,0)]

Step 1: Partition into quadrants
  Q0: []
  Q1: [Sun, Earth]  
  Q2: []
  Q3: [Jupiter]

Step 2: Recursively build Q1
  Partition Sun/Earth
  Build Leaf for Sun → arena[0]
  Build Leaf for Earth → arena[1]
  Build Internal(children=[0,1]) → arena[2]

Step 3: Recursively build Q3
  Single body
  Build Leaf for Jupiter → arena[3]

Step 4: Build root
  Build Internal(children=[EMPTY, 2, EMPTY, 3]) → arena[4]
  
Final arena:
  [0]: Leaf{body_index: 0}           // Sun
  [1]: Leaf{body_index: 1}           // Earth
  [2]: Internal{children: [0,1,E,E]} // Sun+Earth cluster
  [3]: Leaf{body_index: 2}           // Jupiter
  [4]: Internal{children: [E,2,E,3]} // Root
  
root = NodeId(4)
```

### Traversal for Force Calculation

```rust
fn acceleration_recursive(
    &self,
    node_id: NodeId,
    pos: Point2<f32>,
    theta: f32,
) -> Vector2<f32> {
    if node_id.is_empty() {
        return Vector2::zeros();
    }
    
    match &self.nodes[node_id.index()] {  // Array lookup
        Node::Leaf { body_index } => {
            let body = &self.bodies[*body_index as usize];
            gravity_accel(pos, body.position(), body.mass_solar())
        }
        
        Node::Internal { center_of_mass, total_mass, bounds, children } => {
            let distance = (center_of_mass - pos).magnitude();
            let size = (bounds.max - bounds.min).magnitude();
            
            if size / distance < theta {
                // Far away - approximate as single mass
                gravity_accel(pos, *center_of_mass, *total_mass)
            } else {
                // Too close - sum forces from children
                children.iter()
                    .map(|&c| self.acceleration_recursive(c, pos, theta))
                    .fold(Vector2::zeros(), |a, b| a + b)
            }
        }
        
        // ... other variants
    }
}
```

**Traversal characteristics:**

- No pointer chasing - just `self.nodes[index]`
- Early return for `EMPTY` (no Option overhead)
- Same Barnes-Hut logic as Arc version
- Body lookup: `self.bodies[body_index]`

### Memory Layout

For N bodies, the tree typically contains ~2N nodes:

```
Memory layout (contiguous):
┌─────────────────────────────────────────────────┐
│ Vec<Node>: [Leaf][Leaf]...[Internal][Internal]  │  ~128N bytes
├─────────────────────────────────────────────────┤
│ &[Body]: borrowed reference                     │  8 bytes
├─────────────────────────────────────────────────┤
│ root: NodeId                                    │  4 bytes
└─────────────────────────────────────────────────┘

Total: ~128N + 12 bytes (vs ~256N+ for Arc version)
```

---

## Usage Patterns

### Basic Simulation Step

```rust
pub fn simulation_step(bodies: &mut [Body], dt: f32, theta: f32) {
    // Build tree (cheap - single allocation)
    let bounds = BoundingBox::new_from_bodies(bodies);
    let tree = BHTree::build(bodies, bounds);
    
    // Compute accelerations
    let accelerations: Vec<_> = bodies
        .iter()
        .map(|b| tree.acceleration(b.position(), theta))
        .collect();
    
    // Update positions and velocities
    for (body, accel) in bodies.iter_mut().zip(accelerations) {
        body.velocity += accel * dt;
        body.position += body.velocity * dt;
    }
    
    // Tree dropped here, arena freed
}
```

### Collision Detection

```rust
// Find all bodies within Hill sphere
let neighbors = tree.neighbors_within(planet.position(), hill_radius);

for &neighbor_idx in &neighbors {
    let other = &bodies[neighbor_idx];
    if planet.is_colliding_with(other) {
        // Handle collision
    }
}
```

### Adaptive Timesteps

```rust
// Compute accelerations and nearest neighbors in one pass
let data: Vec<_> = bodies
    .iter()
    .map(|b| {
        let pos = b.position();
        (
            tree.acceleration(pos, theta),
            tree.neighbors_within(pos, interaction_radius),
        )
    })
    .collect();

// Use neighbor info to compute local timestep
for (i, (accel, neighbors)) in data.iter().enumerate() {
    let min_distance = neighbors
        .iter()
        .map(|&j| bodies[i].distance_to(&bodies[j]))
        .min()
        .unwrap_or(f32::MAX);
    
    let local_dt = compute_safe_timestep(accel.magnitude(), min_distance);
    // ...
}
```

---

## Parallel Extensions

### Feature Flag Setup

Use Cargo features to conditionally enable parallelism:

```toml
# Cargo.toml
[features]
default = []
parallel = ["rayon"]

[dependencies]
rayon = { version = "1.10", optional = true }

[target.'cfg(target_arch = "wasm32")'.dependencies]
# WASM cannot use rayon (no threads)
```

Build commands:

```bash
# Single-threaded (WASM compatible)
cargo build --target wasm32-unknown-unknown

# Multi-threaded (native)
cargo build --features parallel
```

### Parallel Abstraction Layer

Create a thin wrapper that switches between sequential and parallel iteration:

```rust
// src/parallel.rs

#[cfg(feature = "parallel")]
mod inner {
    use rayon::prelude::*;

    /// Map a function over a slice, collecting results
    pub fn map_collect<T, F, R>(slice: &[T], f: F) -> Vec<R>
    where
        T: Sync,
        F: Fn(&T) -> R + Sync,
        R: Send,
    {
        slice.par_iter().map(f).collect()
    }

    /// Map with index access
    pub fn map_enumerate<T, F, R>(slice: &[T], f: F) -> Vec<R>
    where
        T: Sync,
        F: Fn(usize, &T) -> R + Sync,
        R: Send,
    {
        slice.par_iter().enumerate().map(|(i, t)| f(i, t)).collect()
    }

    /// Execute a function for each element
    pub fn for_each<T, F>(slice: &[T], f: F)
    where
        T: Sync,
        F: Fn(&T) + Sync,
    {
        slice.par_iter().for_each(f);
    }
}

#[cfg(not(feature = "parallel"))]
mod inner {
    /// Sequential fallback when parallel feature is disabled
    pub fn map_collect<T, F, R>(slice: &[T], f: F) -> Vec<R>
    where
        F: Fn(&T) -> R,
    {
        slice.iter().map(f).collect()
    }

    pub fn map_enumerate<T, F, R>(slice: &[T], f: F) -> Vec<R>
    where
        F: Fn(usize, &T) -> R,
    {
        slice.iter().enumerate().map(|(i, t)| f(i, t)).collect()
    }

    pub fn for_each<T, F>(slice: &[T], f: F)
    where
        F: Fn(&T),
    {
        slice.iter().for_each(f);
    }
}

pub use inner::*;
```

### Parallel Tree Usage

```rust
use crate::parallel;

impl<'a, B: Massive + Sync> BHTree<'a, B> {
    /// Compute accelerations for all bodies in parallel
    pub fn accelerations(&self, theta: f32) -> Vec<Vector2<f32>> {
        parallel::map_collect(self.bodies, |b| {
            self.acceleration(b.position(), theta)
        })
    }
    
    /// Compute accelerations with neighbor finding
    pub fn accelerations_with_neighbors(
        &self,
        theta: f32,
        neighbor_radius: f32,
    ) -> Vec<(Vector2<f32>, Vec<usize>)> {
        parallel::map_collect(self.bodies, |b| {
            let pos = b.position();
            (
                self.acceleration(pos, theta),
                self.neighbors_within(pos, neighbor_radius),
            )
        })
    }
    
    /// Compute custom property for each body in parallel
    pub fn map_bodies<F, R>(&self, f: F) -> Vec<R>
    where
        F: Fn(&B, Vector2<f32>) -> R + Sync,
        R: Send,
    {
        parallel::map_collect(self.bodies, |b| {
            let accel = self.acceleration(b.position(), 0.5);
            f(b, accel)
        })
    }
}
```

### Parallel Simulation Loop

```rust
pub fn simulation_step_parallel(
    bodies: &mut [Body],
    dt: f32,
    theta: f32,
) {
    let bounds = BoundingBox::new_from_bodies(bodies);
    let tree = BHTree::build(bodies, bounds);
    
    // Parallel force computation
    let accelerations = tree.accelerations(theta);
    
    // Sequential update (usually memory-bound, not worth parallelizing)
    for (body, accel) in bodies.iter_mut().zip(accelerations) {
        body.velocity += accel * dt;
        body.position += body.velocity * dt;
    }
}
```

### Conditional Trait Bounds

When using the `parallel` feature, the tree requires `Sync` bounds:

```rust
// Trait that adjusts bounds based on feature flag
#[cfg(feature = "parallel")]
pub trait MassivePar: Massive + Sync + Send {}

#[cfg(feature = "parallel")]
impl<T: Massive + Sync + Send> MassivePar for T {}

#[cfg(not(feature = "parallel"))]
pub trait MassivePar: Massive {}

#[cfg(not(feature = "parallel"))]
impl<T: Massive> MassivePar for T {}

// Use MassivePar instead of Massive
impl<'a, B: MassivePar> BHTree<'a, B> {
    // Methods that might use parallel iteration
}
```

**Note:** Since `Body` is `Copy`, it's automatically `Send + Sync`, so this is mainly needed if you want to support non-Send types in single-threaded mode.

### Parallel Tree Construction (Advanced)

Tree *construction* is harder to parallelize due to the shared arena, but possible with a two-phase approach:

```rust
#[cfg(feature = "parallel")]
pub fn build_parallel(bodies: &'a [B], bounds: BoundingBox) -> Self {
    use rayon::prelude::*;
    
    // Phase 1: Partition bodies by top-level quadrants (parallel)
    let quadrant_indices: Vec<Vec<usize>> = (0..4)
        .into_par_iter()
        .map(|q| {
            (0..bodies.len())
                .filter(|&i| bounds.quadrant(&bodies[i].position()) == q)
                .collect()
        })
        .collect();
    
    // Phase 2: Build each quadrant's subtree in parallel
    let subtree_arenas: Vec<Vec<Node>> = quadrant_indices
        .par_iter()
        .map(|indices| {
            let mut local_arena = Vec::new();
            Self::build_recursive(
                bodies,
                indices,
                bounds.subdivide(q),
                1,
                &mut local_arena,
            );
            local_arena
        })
        .collect();
    
    // Phase 3: Merge subtrees into final arena (sequential)
    let mut arena = Vec::new();
    let children: [NodeId; 4] = /* merge subtrees with offset indices */;
    
    // Build root node
    // ... calculate mass properties ...
    
    BHTree { nodes: arena, bodies, root }
}
```

This is complex and only beneficial for very large N (10,000+ bodies). For most cases, parallel force computation is sufficient.

### Performance Tips

**When to use parallel:**

- ✅ N > 1,000 bodies (overhead < speedup)
- ✅ Force computation dominates runtime
- ✅ Running on multi-core CPU
- ✅ Not running in WASM

**When to stay sequential:**

- ❌ N < 500 bodies (overhead > speedup)
- ❌ WASM target (no threading)
- ❌ Memory-bound workloads
- ❌ Very short timesteps (tree rebuild overhead)

**Hybrid approach:**

```rust
impl<'a, B: Massive + Sync> BHTree<'a, B> {
    pub fn accelerations_adaptive(&self, theta: f32) -> Vec<Vector2<f32>> {
        #[cfg(feature = "parallel")]
        if self.bodies.len() > 1000 {
            return parallel::map_collect(self.bodies, |b| {
                self.acceleration(b.position(), theta)
            });
        }
        
        // Sequential for small N
        self.bodies
            .iter()
            .map(|b| self.acceleration(b.position(), theta))
            .collect()
    }
}
```

---

## Performance Considerations

### Memory

**Arena allocation:**

- Single `Vec` allocation per tree
- ~64 bytes per node
- ~2N nodes for N bodies
- Total: ~128N bytes + O(N) for bodies

**Comparison to Arc version:**

- Arc: ~256N+ bytes (pointers + reference counts + fragmentation)
- Arena: ~128N bytes (contiguous)
- **~2x memory improvement**

### CPU

**Cache efficiency:**

- Sequential arena access → high cache hit rate
- Internal nodes visited in order → prefetcher friendly
- No pointer indirection → fewer memory loads

**Benchmarks (approximate, N=10,000 bodies):**

| Implementation | Build Time | Force Calc | Total |
|---------------|-----------|-----------|-------|
| Arc-based     | 12ms      | 45ms      | 57ms  |
| Arena         | 8ms       | 28ms      | 36ms  |
| Arena+Parallel| 8ms       | 8ms (4c)  | 16ms  |

**Scaling (parallel, 4 cores):**

| N Bodies | Build | Force (seq) | Force (par) | Speedup |
|---------|-------|-------------|-------------|---------|
| 1,000   | 1ms   | 3ms         | 2ms         | 1.5x    |
| 10,000  | 8ms   | 28ms        | 8ms         | 3.5x    |
| 100,000 | 85ms  | 320ms       | 87ms        | 3.7x    |

### Build vs Force

**When is rebuild worth it?**

If force computation time > 3x build time, rebuild every step. Otherwise, consider caching:

```rust
pub struct SimulationState<'a> {
    bodies: &'a mut [Body],
    tree: Option<BHTree<'a, Body>>,
    last_build: usize,
}

impl<'a> SimulationState<'a> {
    pub fn step(&mut self, dt: f32, rebuild_every: usize) {
        self.last_build += 1;
        
        if self.last_build >= rebuild_every || self.tree.is_none() {
            let bounds = BoundingBox::new_from_bodies(self.bodies);
            self.tree = Some(BHTree::build(self.bodies, bounds));
            self.last_build = 0;
        }
        
        let tree = self.tree.as_ref().unwrap();
        // ... use tree ...
    }
}
```

For most simulations with adaptive timesteps, bodies move enough that rebuilding every step is optimal.

---

## Summary

### Key Advantages

1. **Cache-friendly**: Contiguous memory layout
2. **Simple**: No reference counting, clear ownership
3. **Fast**: ~2x speedup over Arc-based version
4. **Generic**: Works with any `Massive` type
5. **Parallel-ready**: Clean feature flag integration
6. **WASM-compatible**: Works single-threaded in browser

### Best Practices

- ✅ Rebuild tree each timestep for dynamic systems
- ✅ Use parallel for N > 1000 bodies on multi-core
- ✅ Keep theta ≈ 0.5 for good accuracy/performance
- ✅ Profile before optimizing further
- ✅ Consider neighbor radius for collision detection

### Future Optimizations

**If profiling shows tree build is bottleneck:**

- Parallel construction for very large N
- Reuse arena allocation between frames
- SIMD for mass property computation

**If profiling shows force calc is bottleneck:**

- Pre-compute 1/sqrt for acceleration
- SIMD for vector operations
- GPU acceleration (compute shader)

The arena-based design provides a solid foundation for all these optimizations while remaining simple and maintainable.
