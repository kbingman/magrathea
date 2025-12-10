//! Arena-based Barnes-Hut quadtree for efficient N-body simulations.
//!
//! This module provides a cache-friendly Barnes-Hut implementation using arena allocation.
//! Instead of using `Arc` pointers, nodes are stored contiguously in a `Vec` and reference
//! each other by index, improving cache locality and reducing allocation overhead.
//!
//! # Key Features
//!
//! - **Cache-friendly**: All nodes stored contiguously in memory
//! - **Generic**: Works with any type implementing the `Massive` trait
//! - **Functional build**: Each `build()` creates a fresh arena (no mutation)
//! - **Zero-copy**: References original body slice rather than cloning
//!
//! # Example
//!
//! ```rust
//! use nalgebra::Point2;
//! use nbody::body::Body;
//! use nbody::arena_bhtree::{BHTree, BoundingBox};
//!
//! let bodies = vec![
//!     Body::new_solar_masses(1.0, [0.0, 0.0], [0.0, 0.0]),
//!     Body::new_earth_masses(1.0, [1.0, 0.0], [0.0, 29.78]),
//! ];
//!
//! let bounds = BoundingBox::new_from_bodies(&bodies);
//! let tree = BHTree::build(&bodies, bounds);
//!
//! // Compute acceleration on Earth from all other bodies
//! let earth_pos = Point2::new(1.0, 0.0);
//! let accel = tree.acceleration(earth_pos, 0.5);
//! ```

use nalgebra::{Point2, Vector2};
use units::mass::Mass;

/// Gravitational constant in astronomical units (AU³ solar_mass⁻¹ year⁻²)
/// G = 4π² ≈ 39.478 in these units
const G: f64 = 39.478417;

/// Softening parameter squared to prevent gravitational singularities
const SOFTENING_SQ: f64 = 1e-10;

/// Trait for objects that can participate in Barnes-Hut gravitational calculations.
///
/// This trait provides the minimal interface needed for the Barnes-Hut algorithm:
/// position and mass. Types implementing this trait must be `Copy` to avoid
/// expensive cloning during tree construction.
///
/// # Examples
///
/// ```rust
/// use nalgebra::Point2;
/// use units::mass::Mass;
/// use nbody::arena_bhtree::Massive;
///
/// #[derive(Clone, Copy)]
/// struct Particle {
///     pos: Point2<f64>,
///     mass: Mass,
/// }
///
/// impl Massive for Particle {
///     fn position(&self) -> Point2<f64> {
///         self.pos
///     }
///
///     fn mass(&self) -> Mass {
///         self.mass
///     }
/// }
/// ```
pub trait Massive: Copy {
    /// Returns the position of this body in AU coordinates
    fn position(&self) -> Point2<f64>;

    /// Returns the mass of this body
    fn mass(&self) -> Mass;

    /// Returns mass in solar masses as f64 (cached to avoid repeated conversion)
    fn mass_solar(&self) -> f64 {
        self.mass().to_solar_masses()
    }
}

/// A rectangular bounding box in 2D space using raw AU coordinates.
///
/// Used by the quadtree to define spatial boundaries for recursive subdivision
/// and to determine which quadrant bodies belong to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingBox {
    /// Minimum corner (bottom-left) of the bounding box in AU
    pub min: Point2<f64>,
    /// Maximum corner (top-right) of the bounding box in AU
    pub max: Point2<f64>,
}

impl BoundingBox {
    /// Creates a bounding box that encompasses all the given bodies.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use nbody::body::Body;
    /// use nbody::arena_bhtree::BoundingBox;
    ///
    /// let bodies = vec![
    ///     Body::new_solar_masses(1.0, [-2.0, -1.0], [0.0, 0.0]),
    ///     Body::new_earth_masses(1.0, [3.0, 4.0], [0.0, 0.0]),
    /// ];
    ///
    /// let bounds = BoundingBox::new_from_bodies(&bodies);
    /// ```
    pub fn new_from_bodies<B: Massive>(bodies: &[B]) -> Self {
        bodies.iter().fold(
            Self {
                min: Point2::new(f64::INFINITY, f64::INFINITY),
                max: Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            },
            |bounds, body| {
                let pos = body.position();
                Self {
                    min: Point2::new(bounds.min.x.min(pos.x), bounds.min.y.min(pos.y)),
                    max: Point2::new(bounds.max.x.max(pos.x), bounds.max.y.max(pos.y)),
                }
            },
        )
    }

    /// Returns the center point of the bounding box
    fn center(&self) -> Point2<f64> {
        Point2::new(
            (self.min.x + self.max.x) / 2.0,
            (self.min.y + self.max.y) / 2.0,
        )
    }

    /// Determines which quadrant (0-3) a point belongs to within this bounding box.
    ///
    /// Quadrant layout:
    /// ```text
    /// +-------+-------+
    /// |   2   |   3   |  (top-left, top-right)
    /// +-------+-------+
    /// |   0   |   1   |  (bottom-left, bottom-right)
    /// +-------+-------+
    /// ```
    fn quadrant(&self, point: &Point2<f64>) -> usize {
        let center = self.center();
        let x_bit = (point.x > center.x) as usize;
        let y_bit = (point.y > center.y) as usize;
        x_bit | (y_bit << 1)
    }

    /// Creates a sub-bounding box for the specified quadrant (0-3)
    fn subdivide(&self, quadrant: usize) -> Self {
        let center = self.center();
        let min = Point2::new(
            if quadrant & 1 != 0 {
                center.x
            } else {
                self.min.x
            },
            if quadrant & 2 != 0 {
                center.y
            } else {
                self.min.y
            },
        );
        let max = Point2::new(
            if quadrant & 1 != 0 {
                self.max.x
            } else {
                center.x
            },
            if quadrant & 2 != 0 {
                self.max.y
            } else {
                center.y
            },
        );
        BoundingBox { min, max }
    }
}

/// Index into the node arena.
///
/// Uses `u32` instead of `usize` to save space (4 bytes vs 8 on 64-bit systems).
/// Each internal node stores 4 children, so this saves 16 bytes per node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeId(u32);

impl NodeId {
    /// Sentinel value representing an empty/null node
    pub const EMPTY: NodeId = NodeId(u32::MAX);

    /// Creates a new NodeId from an index
    fn new(index: usize) -> Self {
        debug_assert!(index < u32::MAX as usize, "NodeId overflow");
        NodeId(index as u32)
    }

    /// Returns the index into the arena
    fn index(self) -> usize {
        self.0 as usize
    }

    /// Checks if this is the empty sentinel
    pub(crate) fn is_empty(self) -> bool {
        self == Self::EMPTY
    }
}

/// A node in the arena-based quadtree.
///
/// Nodes reference each other by index (`NodeId`) rather than pointers,
/// allowing all nodes to be stored contiguously for better cache locality.
#[derive(Clone, Copy)]
pub enum Node {
    /// No bodies in this region
    Empty,

    /// Single body (stores index into original body slice)
    Leaf {
        /// Index into the original `&[B]` slice
        body_index: u32,
    },

    /// Multiple bodies at maximum tree depth
    LeafMulti {
        /// Starting index in the original body slice
        start: u32,
        /// Number of bodies in this leaf
        count: u16,
    },

    /// Internal node representing a region of space
    Internal {
        /// Center of mass of all bodies in this subtree
        center_of_mass: Point2<f64>,
        /// Total mass in solar masses (pre-converted for performance)
        total_mass: f64,
        /// Spatial bounds of this node
        bounds: BoundingBox,
        /// Four child quadrants [bottom-left, bottom-right, top-left, top-right]
        children: [NodeId; 4],
    },
}

/// Arena-based Barnes-Hut quadtree.
///
/// All nodes are stored contiguously in a `Vec` and reference each other by index.
/// The tree does not own the bodies - it only stores indices into the original slice.
/// This makes rebuilding cheap and avoids unnecessary copying.
///
/// # Lifetime
///
/// The lifetime `'a` ties the tree to the body slice. The tree is only valid
/// as long as the original bodies exist and remain unchanged.
pub struct BHTree<'a, B: Massive> {
    /// Contiguous storage for all nodes
    nodes: Vec<Node>,
    /// Reference to the original body slice
    bodies: &'a [B],
    /// Root node of the tree
    pub(crate) root: NodeId,
}

impl<'a, B: Massive> BHTree<'a, B> {
    /// Builds a Barnes-Hut quadtree from a collection of bodies.
    ///
    /// This creates a fresh arena and recursively subdivides space until each
    /// leaf contains at most one body (or multiple bodies at maximum depth).
    ///
    /// # Arguments
    ///
    /// * `bodies` - Slice of bodies to organize in the quadtree
    /// * `bounds` - Spatial boundaries for the root node
    ///
    /// # Examples
    ///
    /// ```rust
    /// use nbody::body::Body;
    /// use nbody::arena_bhtree::{BHTree, BoundingBox};
    ///
    /// let bodies = vec![
    ///     Body::new_solar_masses(1.0, [0.0, 0.0], [0.0, 0.0]),
    ///     Body::new_earth_masses(1.0, [1.0, 0.0], [0.0, 29.78]),
    /// ];
    ///
    /// let bounds = BoundingBox::new_from_bodies(&bodies);
    /// let tree = BHTree::build(&bodies, bounds);
    /// ```
    pub fn build(bodies: &'a [B], bounds: BoundingBox) -> Self {
        let mut arena = Vec::with_capacity(bodies.len() * 2);
        let indices: Vec<usize> = (0..bodies.len()).collect();

        let root = Self::build_recursive(bodies, &indices, bounds, 0, &mut arena);

        BHTree {
            nodes: arena,
            bodies,
            root,
        }
    }

    /// Recursively builds the quadtree with depth tracking to prevent infinite recursion
    fn build_recursive(
        bodies: &[B],
        indices: &[usize],
        bounds: BoundingBox,
        depth: usize,
        arena: &mut Vec<Node>,
    ) -> NodeId {
        const MAX_DEPTH: usize = 30;

        match indices {
            // No bodies - return empty sentinel
            [] => NodeId::EMPTY,

            // Single body - create leaf
            [single] => {
                let id = NodeId::new(arena.len());
                arena.push(Node::Leaf {
                    body_index: *single as u32,
                });
                id
            }

            // Max depth reached - store all bodies in a multi-body leaf
            _ if depth >= MAX_DEPTH => {
                let id = NodeId::new(arena.len());
                arena.push(Node::LeafMulti {
                    start: indices[0] as u32,
                    count: indices.len() as u16,
                });
                id
            }

            // Multiple bodies - subdivide into quadrants
            indices => {
                // Partition indices by quadrant
                let mut quadrants: [Vec<usize>; 4] = Default::default();
                for &i in indices {
                    let q = bounds.quadrant(&bodies[i].position());
                    quadrants[q].push(i);
                }

                // Recursively build children
                let children: [NodeId; 4] = std::array::from_fn(|q| {
                    Self::build_recursive(
                        bodies,
                        &quadrants[q],
                        bounds.subdivide(q),
                        depth + 1,
                        arena,
                    )
                });

                // Compute mass properties in a single pass
                let (total_mass, weighted_pos) =
                    indices
                        .iter()
                        .fold((0.0f64, Vector2::zeros()), |(mass, pos), &i| {
                            let b = &bodies[i];
                            let m = b.mass_solar();
                            (mass + m, pos + b.position().coords * m)
                        });

                let center_of_mass = Point2::from(weighted_pos / total_mass);

                // Push internal node and return its ID
                let id = NodeId::new(arena.len());
                arena.push(Node::Internal {
                    center_of_mass,
                    total_mass,
                    bounds,
                    children,
                });
                id
            }
        }
    }

    /// Computes the gravitational acceleration at a point using the Barnes-Hut algorithm.
    ///
    /// This traverses the quadtree and decides whether to approximate distant clusters
    /// or recurse for nearby regions based on the theta parameter.
    ///
    /// # Arguments
    ///
    /// * `pos` - Position to compute acceleration at (AU coordinates)
    /// * `theta` - Approximation parameter (typically 0.5-2.0)
    ///   - Smaller θ = more accurate, slower
    ///   - Larger θ = less accurate, faster
    ///
    /// # Returns
    ///
    /// Gravitational acceleration vector (AU/year²)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use nalgebra::Point2;
    /// use nbody::body::Body;
    /// use nbody::arena_bhtree::{BHTree, BoundingBox};
    ///
    /// let bodies = vec![
    ///     Body::new_solar_masses(1.0, [0.0, 0.0], [0.0, 0.0]),
    ///     Body::new_earth_masses(1.0, [5.0, 0.0], [0.0, 0.0]),
    /// ];
    ///
    /// let bounds = BoundingBox::new_from_bodies(&bodies);
    /// let tree = BHTree::build(&bodies, bounds);
    ///
    /// let pos = Point2::new(5.0, 0.0);
    /// let accel = tree.acceleration(pos, 0.5);
    /// ```
    pub fn acceleration(&self, pos: Point2<f64>, theta: f64) -> Vector2<f64> {
        self.acceleration_recursive(self.root, pos, theta)
    }

    /// Recursive helper for acceleration computation
    fn acceleration_recursive(
        &self,
        node_id: NodeId,
        pos: Point2<f64>,
        theta: f64,
    ) -> Vector2<f64> {
        // Early return for empty nodes
        if node_id.is_empty() {
            return Vector2::zeros();
        }

        match &self.nodes[node_id.index()] {
            Node::Empty => Vector2::zeros(),

            Node::Leaf { body_index } => {
                let body = &self.bodies[*body_index as usize];
                gravity_accel(pos, body.position(), body.mass_solar())
            }

            Node::LeafMulti { start, count } => {
                // Sum accelerations from all bodies in multi-body leaf
                (*start..(*start + *count as u32))
                    .map(|i| {
                        let body = &self.bodies[i as usize];
                        gravity_accel(pos, body.position(), body.mass_solar())
                    })
                    .fold(Vector2::zeros(), |a, b| a + b)
            }

            Node::Internal {
                center_of_mass,
                total_mass,
                bounds,
                children,
            } => {
                // Barnes-Hut criterion: s/d < θ
                let diff = *center_of_mass - pos;
                let distance = diff.magnitude();
                let size = (bounds.max - bounds.min).magnitude();

                if size / distance < theta {
                    // Far enough - use center of mass approximation
                    gravity_accel(pos, *center_of_mass, *total_mass)
                } else {
                    // Too close - recurse into children
                    children
                        .iter()
                        .map(|&child| self.acceleration_recursive(child, pos, theta))
                        .fold(Vector2::zeros(), |a, b| a + b)
                }
            }
        }
    }

    /// Finds all bodies within a given radius of a target position.
    ///
    /// This performs a proximity query on the quadtree, pruning entire subtrees
    /// whose bounding boxes don't intersect the search circle.
    ///
    /// # Arguments
    ///
    /// * `pos` - Center position for the search (AU coordinates)
    /// * `radius` - Search radius (AU)
    ///
    /// # Returns
    ///
    /// Vector of body indices within the search radius
    ///
    /// # Examples
    ///
    /// ```rust
    /// use nalgebra::Point2;
    /// use nbody::body::Body;
    /// use nbody::arena_bhtree::{BHTree, BoundingBox};
    ///
    /// let bodies = vec![
    ///     Body::new_solar_masses(1.0, [0.0, 0.0], [0.0, 0.0]),
    ///     Body::new_earth_masses(1.0, [1.0, 0.0], [0.0, 0.0]),
    ///     Body::new_earth_masses(1.0, [5.0, 0.0], [0.0, 0.0]),
    /// ];
    ///
    /// let bounds = BoundingBox::new_from_bodies(&bodies);
    /// let tree = BHTree::build(&bodies, bounds);
    ///
    /// let neighbors = tree.neighbors_within(Point2::new(0.0, 0.0), 2.0);
    /// assert_eq!(neighbors.len(), 2); // Sun and Earth, not the distant body
    /// ```
    pub fn neighbors_within(&self, pos: Point2<f64>, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        let radius_sq = radius * radius;
        self.neighbors_recursive(self.root, pos, radius_sq, &mut result);
        result
    }

    /// Recursive helper for neighbor finding
    fn neighbors_recursive(
        &self,
        node_id: NodeId,
        pos: Point2<f64>,
        radius_sq: f64,
        result: &mut Vec<usize>,
    ) {
        if node_id.is_empty() {
            return;
        }

        match &self.nodes[node_id.index()] {
            Node::Empty => {}

            Node::Leaf { body_index } => {
                let body = &self.bodies[*body_index as usize];
                if (body.position() - pos).magnitude_squared() <= radius_sq {
                    result.push(*body_index as usize);
                }
            }

            Node::LeafMulti { start, count } => {
                for i in *start..(*start + *count as u32) {
                    let body = &self.bodies[i as usize];
                    if (body.position() - pos).magnitude_squared() <= radius_sq {
                        result.push(i as usize);
                    }
                }
            }

            Node::Internal {
                bounds, children, ..
            } => {
                // Check if bounding box intersects search circle
                if box_intersects_circle(bounds, pos, radius_sq.sqrt()) {
                    for &child in children {
                        self.neighbors_recursive(child, pos, radius_sq, result);
                    }
                }
            }
        }
    }

    /// Returns the number of nodes in the tree (for diagnostics)
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the root node ID (for advanced use cases)
    pub fn root(&self) -> NodeId {
        self.root
    }
}

/// Computes gravitational acceleration toward a mass point.
///
/// Uses Newton's law of universal gravitation with softening to prevent singularities.
///
/// # Arguments
///
/// * `from` - Position experiencing the acceleration
/// * `toward` - Position of the mass
/// * `mass_solar` - Mass in solar masses
///
/// # Returns
///
/// Acceleration vector in AU/year²
#[inline]
fn gravity_accel(from: Point2<f64>, toward: Point2<f64>, mass_solar: f64) -> Vector2<f64> {
    let diff = toward - from;
    let dist_sq = diff.magnitude_squared() + SOFTENING_SQ;
    let dist = dist_sq.sqrt();

    // a = G * M * r / r³ = G * M * r / (r²)^(3/2)
    diff * (G * mass_solar / (dist_sq * dist))
}

/// Checks if a bounding box intersects with a circle.
///
/// Used to prune quadtree branches during proximity queries.
///
/// # Arguments
///
/// * `bounds` - The rectangular bounding box
/// * `center` - Center of the search circle
/// * `radius` - Radius of the search circle
///
/// # Returns
///
/// `true` if the bounding box intersects or is contained within the circle
fn box_intersects_circle(bounds: &BoundingBox, center: Point2<f64>, radius: f64) -> bool {
    // Find closest point on bounding box to circle center
    let closest_x = center.x.clamp(bounds.min.x, bounds.max.x);
    let closest_y = center.y.clamp(bounds.min.y, bounds.max.y);

    // Calculate squared distance from circle center to closest point
    let dx = center.x - closest_x;
    let dy = center.y - closest_y;
    let dist_sq = dx * dx + dy * dy;

    // Circle intersects box if closest point is within radius
    dist_sq <= radius * radius
}
