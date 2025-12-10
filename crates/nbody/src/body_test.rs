use nalgebra::{Point2, Vector2};
use units::{Length, Mass};

use crate::body::{Body, BodyId};

#[test]
fn test_new_solar_masses() {
    let body = Body::new_solar_masses(1.0, [0.0, 0.0], [0.0, 0.0]);

    assert_eq!(body.mass, 1.0);
    assert_eq!(body.position, Point2::new(0.0, 0.0));
    assert_eq!(body.velocity, Vector2::new(0.0, 0.0));
    assert_eq!(body.radius, Length::from_solar_radii(1.0).to_au());
}

#[test]
fn test_new_earth_masses() {
    let body = Body::new_earth_masses(1.0, [1.0, 0.0], [0.0, 6.28]);

    let expected_mass = Mass::from_earth_masses(1.0).to_solar_masses();
    assert!((body.mass - expected_mass).abs() < 1e-10);
    assert_eq!(body.position, Point2::new(1.0, 0.0));
    assert_eq!(body.velocity, Vector2::new(0.0, 6.28));
    assert_eq!(body.radius, Length::from_earth_radii(1.0).to_au());
}

#[test]
fn test_momentum() {
    let body = Body {
        id: BodyId(0),
        mass: 2.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(3.0, 4.0),
    };

    let momentum = body.momentum();
    assert_eq!(momentum, Vector2::new(6.0, 8.0));
}

#[test]
fn test_kinetic_energy() {
    let body = Body {
        id: BodyId(0),
        mass: 2.0,
        radius: 0.01,
        position: Point2::new(0.0, 0.0),
        velocity: Vector2::new(3.0, 4.0),
    };

    // KE = 0.5 * m * v²
    // v² = 3² + 4² = 25
    // KE = 0.5 * 2 * 25 = 25
    let ke = body.kinetic_energy();
    assert_eq!(ke, 25.0);
}

#[test]
fn test_distance_to() {
    let body_a = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(0.0, 0.0),
        velocity: Vector2::new(0.0, 0.0),
    };

    let body_b = Body {
        id: BodyId(1),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(3.0, 4.0),
        velocity: Vector2::new(0.0, 0.0),
    };

    // Distance = sqrt(3² + 4²) = 5
    let distance = body_a.distance_to(&body_b);
    assert_eq!(distance, 5.0);
}

#[test]
fn test_orbital_radius() {
    let body = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(3.0, 4.0),
        velocity: Vector2::new(0.0, 0.0),
    };

    // r = sqrt(3² + 4²) = 5
    let r = body.orbital_radius();
    assert_eq!(r, 5.0);
}

#[test]
fn test_orbital_velocity() {
    let body = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(0.0, 0.0),
        velocity: Vector2::new(3.0, 4.0),
    };

    // v = sqrt(3² + 4²) = 5
    let v = body.orbital_velocity();
    assert_eq!(v, 5.0);
}

#[test]
fn test_specific_angular_momentum() {
    let body = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(0.0, 2.0),
    };

    // L_z = x * v_y - y * v_x = 1.0 * 2.0 - 0.0 * 0.0 = 2.0
    let l = body.specific_angular_momentum();
    assert_eq!(l, 2.0);
}

#[test]
fn test_specific_angular_momentum_zero() {
    let body = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(1.0, 0.0),
    };

    // Radial velocity -> zero angular momentum
    let l = body.specific_angular_momentum();
    assert_eq!(l, 0.0);
}

#[test]
fn test_specific_angular_momentum_negative() {
    let body = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(0.0, 1.0),
        velocity: Vector2::new(2.0, 0.0),
    };

    // L_z = x * v_y - y * v_x = 0.0 * 0.0 - 1.0 * 2.0 = -2.0
    // Negative means clockwise rotation
    let l = body.specific_angular_momentum();
    assert_eq!(l, -2.0);
}

#[test]
fn test_circular_orbit_angular_momentum() {
    // Set up a circular orbit at 1 AU
    const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²
    let mu = G * 1.0; // Sun mass = 1 M☉
    let r = 1.0; // 1 AU
    let v_circular = (mu / r).sqrt();

    let body = Body {
        id: BodyId(0),
        mass: Mass::from_earth_masses(1.0).to_solar_masses(),
        radius: Length::from_earth_radii(1.0).to_au(),
        position: Point2::new(r, 0.0),
        velocity: Vector2::new(0.0, v_circular),
    };

    // For circular orbit: L = r * v
    let expected_l = r * v_circular;
    let actual_l = body.specific_angular_momentum();

    assert!((actual_l - expected_l).abs() / expected_l < 1e-10);
}

#[test]
fn test_massive_trait_position() {
    use crate::arena_bhtree::Massive;

    let body = Body::new_solar_masses(1.0, [5.0, 3.0], [0.0, 0.0]);

    assert_eq!(body.position(), Point2::new(5.0, 3.0));
}

#[test]
fn test_massive_trait_mass() {
    use crate::arena_bhtree::Massive;

    let body = Body::new_earth_masses(1.0, [0.0, 0.0], [0.0, 0.0]);

    let mass = body.mass();
    let expected_mass = Mass::from_earth_masses(1.0);

    // Compare as solar masses
    assert!((mass.to_solar_masses() - expected_mass.to_solar_masses()).abs() < 1e-10);
}

#[test]
fn test_massive_trait_mass_solar() {
    use crate::arena_bhtree::Massive;

    let body = Body::new_solar_masses(2.5, [0.0, 0.0], [0.0, 0.0]);

    assert_eq!(body.mass_solar(), 2.5);
}

#[test]
fn test_body_id_equality() {
    let id1 = BodyId(42);
    let id2 = BodyId(42);
    let id3 = BodyId(43);

    assert_eq!(id1, id2);
    assert_ne!(id1, id3);
}

#[test]
fn test_body_copy() {
    let body1 = Body::new_earth_masses(1.0, [1.0, 2.0], [3.0, 4.0]);
    let body2 = body1; // Should copy, not move

    // Both should be usable
    assert_eq!(body1.mass, body2.mass);
    assert_eq!(body1.position, body2.position);
    assert_eq!(body1.velocity, body2.velocity);
}
