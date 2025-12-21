use super::*;
use nalgebra::{Point2, Vector2};
use units::{Length, Mass};

const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

#[test]
fn circular_orbit_recovers_correct_elements() {
    // Circular orbit at 1 AU around 1 solar mass star
    let pos = Point2::new(1.0, 0.0);
    let vel = Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()); // Circular velocity

    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    assert!(
        (elements.semi_major_axis.to_au() - 1.0).abs() < 0.01,
        "Semi-major axis should be 1.0 AU, got {}",
        elements.semi_major_axis.to_au()
    );
    assert!(
        elements.eccentricity < 0.01,
        "Eccentricity should be near 0, got {}",
        elements.eccentricity
    );
}

#[test]
fn eccentric_orbit_at_periapsis() {
    let a = 1.5;
    let e = 0.3;
    let r_peri = a * (1.0 - e);

    let mu = G * 1.0;
    let v_peri = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();

    let pos = Point2::new(r_peri, 0.0);
    let vel = Vector2::new(0.0, v_peri);

    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    assert!(
        (elements.semi_major_axis.to_au() - a).abs() < 0.01,
        "Semi-major axis should be {} AU, got {}",
        a,
        elements.semi_major_axis.to_au()
    );
    assert!(
        (elements.eccentricity - e).abs() < 0.01,
        "Eccentricity should be {}, got {}",
        e,
        elements.eccentricity
    );
}

#[test]
fn eccentric_orbit_at_apoapsis() {
    let a = 2.0;
    let e = 0.5;
    let r_apo = a * (1.0 + e);

    let mu = G * 1.0;
    let v_apo = (mu * (1.0 - e) / (a * (1.0 + e))).sqrt();

    let pos = Point2::new(r_apo, 0.0);
    let vel = Vector2::new(0.0, v_apo);

    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    assert!(
        (elements.semi_major_axis.to_au() - a).abs() < 0.01,
        "Semi-major axis should be {} AU, got {}",
        a,
        elements.semi_major_axis.to_au()
    );
    assert!(
        (elements.eccentricity - e).abs() < 0.01,
        "Eccentricity should be {}, got {}",
        e,
        elements.eccentricity
    );
}

#[test]
fn arbitrary_orbital_position() {
    // Test at an arbitrary true anomaly
    let a = 1.0;
    let e = 0.2;
    let nu = std::f64::consts::PI / 4.0; // 45 degrees

    let r = a * (1.0 - e * e) / (1.0 + e * nu.cos());
    let mu = G * 1.0;
    let h = (mu * a * (1.0 - e * e)).sqrt(); // Angular momentum

    let pos_x = r * nu.cos();
    let pos_y = r * nu.sin();
    let vel_x = -mu * nu.sin() / h;
    let vel_y = mu * (e + nu.cos()) / h;

    let pos = Point2::new(pos_x, pos_y);
    let vel = Vector2::new(vel_x, vel_y);

    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    assert!(
        (elements.semi_major_axis.to_au() - a).abs() < 0.01,
        "Semi-major axis should be {} AU, got {}",
        a,
        elements.semi_major_axis.to_au()
    );
    assert!(
        (elements.eccentricity - e).abs() < 0.01,
        "Eccentricity should be {}, got {}",
        e,
        elements.eccentricity
    );
}

#[test]
fn period_from_elements() {
    let elements = OrbitalElements {
        semi_major_axis: Length::from_au(1.0),
        eccentricity: 0.0,
        inclination: 0.0,
        longitude_ascending_node: 0.0,
        argument_of_periapsis: 0.0,
        mean_anomaly: 0.0,
    };

    let period = elements.period(Mass::from_solar_masses(1.0));

    // For a = 1 AU, M = 1 M☉: T = 1 year
    assert!(
        (period - 1.0).abs() < 0.01,
        "Period should be 1 year, got {}",
        period
    );
}

#[test]
fn mean_motion_from_elements() {
    let elements = OrbitalElements {
        semi_major_axis: Length::from_au(5.2), // Jupiter's orbit
        eccentricity: 0.0,
        inclination: 0.0,
        longitude_ascending_node: 0.0,
        argument_of_periapsis: 0.0,
        mean_anomaly: 0.0,
    };

    let n = elements.mean_motion(Mass::from_solar_masses(1.0));
    let period = elements.period(Mass::from_solar_masses(1.0));

    // n × T = 2π
    let product = n * period;
    assert!(
        (product - std::f64::consts::TAU).abs() < 0.01,
        "n × T should be 2π, got {}",
        product
    );
}

#[test]
fn round_trip_circular_orbit() {
    // Start with Cartesian
    let pos = Point2::new(1.0, 0.0);
    let vel = Vector2::new(0.0, (G * 1.0 / 1.0).sqrt());

    // Convert to elements
    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    // Convert back to Cartesian
    let (pos2, vel2) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));

    // Should match original
    assert!((pos.x - pos2.x).abs() < 1e-6);
    assert!((pos.y - pos2.y).abs() < 1e-6);
    assert!((vel.x - vel2.x).abs() < 1e-6);
    assert!((vel.y - vel2.y).abs() < 1e-6);
}

#[test]
fn round_trip_eccentric_orbit() {
    // Eccentric orbit at periapsis
    let a = 1.5;
    let e = 0.3;
    let r_peri = a * (1.0 - e);
    let mu = G * 1.0;
    let v_peri = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();

    let pos = Point2::new(r_peri, 0.0);
    let vel = Vector2::new(0.0, v_peri);

    // Convert to elements and back
    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));
    let (pos2, vel2) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));

    // Should match original
    assert!(
        (pos.x - pos2.x).abs() < 1e-6,
        "Position x mismatch: {} vs {}",
        pos.x,
        pos2.x
    );
    assert!(
        (pos.y - pos2.y).abs() < 1e-6,
        "Position y mismatch: {} vs {}",
        pos.y,
        pos2.y
    );
    assert!(
        (vel.x - vel2.x).abs() < 1e-6,
        "Velocity x mismatch: {} vs {}",
        vel.x,
        vel2.x
    );
    assert!(
        (vel.y - vel2.y).abs() < 1e-6,
        "Velocity y mismatch: {} vs {}",
        vel.y,
        vel2.y
    );
}

#[test]
fn round_trip_arbitrary_position() {
    // Random position in orbit
    let a = 2.0;
    let e = 0.5;
    let nu = std::f64::consts::PI / 3.0; // 60 degrees true anomaly

    let r = a * (1.0 - e * e) / (1.0 + e * nu.cos());
    let mu = G * 1.0;
    let h = (mu * a * (1.0 - e * e)).sqrt();

    let pos_x = r * nu.cos();
    let pos_y = r * nu.sin();
    let vel_x = -mu * nu.sin() / h;
    let vel_y = mu * (e + nu.cos()) / h;

    let pos = Point2::new(pos_x, pos_y);
    let vel = Vector2::new(vel_x, vel_y);

    // Round trip
    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));
    let (pos2, vel2) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));

    // Should match
    assert!(
        (pos.x - pos2.x).abs() < 1e-5,
        "Position x: {} vs {}",
        pos.x,
        pos2.x
    );
    assert!(
        (pos.y - pos2.y).abs() < 1e-5,
        "Position y: {} vs {}",
        pos.y,
        pos2.y
    );
    assert!(
        (vel.x - vel2.x).abs() < 1e-5,
        "Velocity x: {} vs {}",
        vel.x,
        vel2.x
    );
    assert!(
        (vel.y - vel2.y).abs() < 1e-5,
        "Velocity y: {} vs {}",
        vel.y,
        vel2.y
    );
}

#[test]
fn round_trip_with_rotation() {
    // Orbit rotated by argument of periapsis
    let elements = OrbitalElements {
        semi_major_axis: Length::from_au(1.0),
        eccentricity: 0.3,
        inclination: 0.0,
        longitude_ascending_node: 0.0,
        argument_of_periapsis: std::f64::consts::PI / 4.0, // 45 degrees
        mean_anomaly: std::f64::consts::PI / 2.0,          // 90 degrees
    };

    // Convert to Cartesian
    let (pos, vel) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));

    // Convert back
    let elements2 = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

    // Check orbital parameters match
    assert!(
        (elements.semi_major_axis.to_au() - elements2.semi_major_axis.to_au()).abs() < 1e-6,
        "Semi-major axis mismatch"
    );
    assert!(
        (elements.eccentricity - elements2.eccentricity).abs() < 1e-6,
        "Eccentricity mismatch"
    );
    // Argument of periapsis might differ by 2π, so check modulo
    let omega_diff = (elements.argument_of_periapsis - elements2.argument_of_periapsis)
        .abs()
        .rem_euclid(std::f64::consts::TAU);
    assert!(
        omega_diff < 1e-6 || (std::f64::consts::TAU - omega_diff).abs() < 1e-6,
        "Argument of periapsis mismatch: {} vs {}",
        elements.argument_of_periapsis,
        elements2.argument_of_periapsis
    );
}

#[test]
fn solve_keplers_equation_at_various_eccentricities() {
    use crate::bodies::orbital_elements::solve_keplers_equation;

    // Test for various eccentricities
    for e in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9] {
        for m in [
            0.0,
            std::f64::consts::PI / 4.0,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI,
        ] {
            let ecc_anom = solve_keplers_equation(m, e);

            // Verify Kepler's equation: M = E - e sin(E)
            let m_check = ecc_anom - e * ecc_anom.sin();

            assert!(
                (m - m_check).abs() < 1e-10,
                "Kepler equation not satisfied for e={}, M={}: got {}",
                e,
                m,
                m_check
            );
        }
    }
}

#[test]
fn elements_to_cartesian_at_periapsis() {
    // Should be at periapsis when mean_anomaly = 0
    let elements = OrbitalElements {
        semi_major_axis: Length::from_au(1.5),
        eccentricity: 0.3,
        inclination: 0.0,
        longitude_ascending_node: 0.0,
        argument_of_periapsis: 0.0,
        mean_anomaly: 0.0,
    };

    let (pos, _vel) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));

    let r = pos.coords.magnitude();
    let r_peri = 1.5 * (1.0 - 0.3);

    assert!(
        (r - r_peri).abs() < 1e-6,
        "Should be at periapsis: r={}, r_peri={}",
        r,
        r_peri
    );
}

#[test]
fn different_stellar_masses() {
    // Test that orbital elements scale correctly with stellar mass
    let pos = Point2::new(1.0, 0.0);

    // For a 2 solar mass star, circular velocity is sqrt(2) times larger
    let vel = Vector2::new(0.0, (G * 2.0 / 1.0).sqrt());

    let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(2.0));

    assert!(
        (elements.semi_major_axis.to_au() - 1.0).abs() < 0.01,
        "Semi-major axis should still be 1.0 AU"
    );
    assert!(
        elements.eccentricity < 0.01,
        "Should still be circular orbit"
    );

    // Period should be 1/sqrt(2) years for 2x mass
    let period = elements.period(Mass::from_solar_masses(2.0));
    let expected_period = 1.0 / 2.0_f64.sqrt();
    assert!(
        (period - expected_period).abs() < 0.01,
        "Period should be {}, got {}",
        expected_period,
        period
    );
}
