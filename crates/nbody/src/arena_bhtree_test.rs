use nalgebra::Point2;
use units::Mass;

use crate::arena_bhtree::{BHTree, BoundingBox, Massive};

#[derive(Clone, Copy)]
struct TestBody {
    pos: Point2<f64>,
    mass: f64,
}

impl Massive for TestBody {
    fn position(&self) -> Point2<f64> {
        self.pos
    }

    fn mass(&self) -> Mass {
        Mass::from_solar_masses(self.mass)
    }

    fn mass_solar(&self) -> f64 {
        self.mass
    }
}

#[test]
fn test_empty_tree() {
    let bodies: Vec<TestBody> = vec![];
    let bounds = BoundingBox {
        min: Point2::new(-10.0, -10.0),
        max: Point2::new(10.0, 10.0),
    };
    let tree = BHTree::build(&bodies, bounds);

    assert_eq!(tree.node_count(), 0);
    assert!(tree.root.is_empty());
}

#[test]
fn test_single_body() {
    let bodies = vec![TestBody {
        pos: Point2::new(1.0, 2.0),
        mass: 1.0,
    }];
    let bounds = BoundingBox::new_from_bodies(&bodies);
    let tree = BHTree::build(&bodies, bounds);

    assert_eq!(tree.node_count(), 1);
}

#[test]
fn test_acceleration_simple() {
    let bodies = vec![
        TestBody {
            pos: Point2::new(0.0, 0.0),
            mass: 1.0,
        },
        TestBody {
            pos: Point2::new(1.0, 0.0),
            mass: 0.000003, // ~Earth mass
        },
    ];

    let bounds = BoundingBox::new_from_bodies(&bodies);
    let tree = BHTree::build(&bodies, bounds);

    let accel = tree.acceleration(Point2::new(1.0, 0.0), 0.5);

    // Earth should be pulled toward Sun (negative x direction)
    assert!(accel.x < 0.0);
    assert!(accel.y.abs() < 1e-6);
}

#[test]
fn test_neighbors() {
    let bodies = vec![
        TestBody {
            pos: Point2::new(0.0, 0.0),
            mass: 1.0,
        },
        TestBody {
            pos: Point2::new(1.0, 0.0),
            mass: 1.0,
        },
        TestBody {
            pos: Point2::new(10.0, 0.0),
            mass: 1.0,
        },
    ];

    let bounds = BoundingBox::new_from_bodies(&bodies);
    let tree = BHTree::build(&bodies, bounds);

    let neighbors = tree.neighbors_within(Point2::new(0.0, 0.0), 2.0);
    assert_eq!(neighbors.len(), 2); // Should find bodies 0 and 1, not 2
}
