use super::*;
use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use units::{Mass, Velocity};

const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

#[test]
fn creates_body_with_correct_orbital_elements() {
    let pos = Point2::new(1.0, 0.0);
    let vel = Vector2::new(0.0, (G * 1.0 / 1.0).sqrt());

    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        pos,
        vel,
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    assert!((body.semi_major_axis.to_au() - 1.0).abs() < 0.01);
    assert!(body.eccentricity < 0.01);
}

#[test]
fn total_mass_is_core_plus_envelope() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(10.0),
        Mass::from_earth_masses(5.0),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::mini_neptune(),
    );

    assert!(
        (body.total_mass().to_earth_masses() - 15.0).abs() < 1e-6,
        "Total mass should be core + envelope"
    );
}

#[test]
fn hill_radius_scales_correctly() {
    // Earth at 1 AU has Hill radius ~0.01 AU
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let r_hill = body.hill_radius(Mass::from_solar_masses(1.0));

    // R_H = a × (M / 3M_star)^(1/3)
    let expected = 1.0 * (Mass::from_earth_masses(1.0).to_solar_masses() / 3.0).powf(1.0 / 3.0);

    assert!(
        (r_hill.to_au() - expected).abs() < 0.001,
        "Hill radius should be ~{} AU, got {}",
        expected,
        r_hill.to_au()
    );
}

#[test]
fn hill_radius_increases_with_mass() {
    let small_body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let large_body = DiscreteBody::new(
        Mass::from_earth_masses(10.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let r_hill_small = small_body.hill_radius(Mass::from_solar_masses(1.0));
    let r_hill_large = large_body.hill_radius(Mass::from_solar_masses(1.0));

    assert!(
        r_hill_large > r_hill_small,
        "More massive body should have larger Hill radius"
    );

    // Should scale as M^(1/3)
    let mass_ratio: f64 = 10.0;
    let expected_ratio = mass_ratio.powf(1.0 / 3.0);
    let actual_ratio = r_hill_large.to_au() / r_hill_small.to_au();

    assert!(
        (actual_ratio - expected_ratio).abs() < 0.01,
        "Hill radius should scale as M^(1/3)"
    );
}

#[test]
fn escape_velocity_for_earth() {
    // Create an Earth-like body
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, 6.28),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let v_esc = body.escape_velocity();

    // Earth's escape velocity is ~11.2 km/s ≈ 2.37 AU/year
    // But our body has estimated radius, so just check it's positive and reasonable
    assert!(
        v_esc.to_au_per_year() > 0.0,
        "Escape velocity should be positive"
    );
    assert!(
        v_esc.to_au_per_year() < 10.0,
        "Escape velocity should be reasonable (<10 AU/year)"
    );
}

#[test]
fn gravitational_focusing_enhances_accretion() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    // Low velocity dispersion → strong focusing
    let low_sigma = Velocity::from_meters_per_sec(10.0);
    let f_g_low = body.gravitational_focusing_factor(low_sigma);

    // High velocity dispersion → weak focusing
    let high_sigma = Velocity::from_meters_per_sec(1000.0);
    let f_g_high = body.gravitational_focusing_factor(high_sigma);

    assert!(
        f_g_low > f_g_high,
        "Lower velocity dispersion should give stronger focusing"
    );
    assert!(f_g_low >= 1.0, "Focusing factor must be ≥ 1");
    assert!(f_g_high >= 1.0, "Focusing factor must be ≥ 1");
}

#[test]
fn focusing_factor_formula() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(0.1),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma = Velocity::from_au_per_year(1.0);
    let f_g = body.gravitational_focusing_factor(sigma);

    // F_g = 1 + (v_esc / σ)²
    let v_esc = body.escape_velocity().to_au_per_year();
    let expected = 1.0 + (v_esc / sigma.to_au_per_year()).powi(2);

    assert!(
        (f_g - expected).abs() < 1e-6,
        "Focusing factor should match formula: expected {}, got {}",
        expected,
        f_g
    );
}

#[test]
fn orbital_frequency_matches_kepler() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let omega = body.orbital_frequency(Mass::from_solar_masses(1.0));
    let period = body.orbital_period(Mass::from_solar_masses(1.0));

    // ω × T = 2π
    let product = omega * period;
    assert!(
        (product - std::f64::consts::TAU).abs() < 0.01,
        "ω × T should equal 2π"
    );
}

#[test]
fn update_orbital_elements_after_position_change() {
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    // Initially at a = 1 AU
    assert!((body.semi_major_axis.to_au() - 1.0).abs() < 0.01);

    // Move to different orbit
    body.position = Point2::new(5.0, 0.0);
    body.velocity = Vector2::new(0.0, (G * 1.0 / 5.0).sqrt());

    // Update elements
    body.update_orbital_elements(Mass::from_solar_masses(1.0));

    // Should now be at a = 5 AU
    assert!(
        (body.semi_major_axis.to_au() - 5.0).abs() < 0.01,
        "Semi-major axis should update to new orbit"
    );
}

#[test]
fn conversion_to_nbody_body() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(10.0),
        Mass::from_earth_masses(5.0),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::mini_neptune(),
    );

    let nbody_body = body.to_nbody_body();

    assert_eq!(nbody_body.position, body.position);
    assert_eq!(nbody_body.velocity, body.velocity);
    assert_eq!(
        nbody_body.mass,
        body.total_mass().to_solar_masses(),
        "nbody body should use total mass"
    );
}

#[test]
fn as_massive_creates_lightweight_body() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::from_earth_masses(0.5),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, 6.28),
        Mass::from_solar_masses(1.0),
        Composition::mini_neptune(),
    );

    // Create lightweight MassiveBody for tree calculations
    let massive = body.as_massive();

    use nbody::arena_bhtree::Massive;

    assert_eq!(massive.position(), body.position);
    assert_eq!(massive.mass(), body.total_mass());
    assert_eq!(massive.mass_solar(), body.total_mass().to_solar_masses());
}

#[test]
fn envelope_state_none_for_no_envelope() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, 6.28),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    assert!(
        matches!(body.envelope_state, EnvelopeState::None),
        "Body with no envelope mass should have EnvelopeState::None"
    );
}

#[test]
fn envelope_state_hydrostatic_with_envelope() {
    let envelope_mass = Mass::from_earth_masses(2.0);
    let body = DiscreteBody::new(
        Mass::from_earth_masses(10.0),
        envelope_mass,
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::mini_neptune(),
    );

    match body.envelope_state {
        EnvelopeState::Hydrostatic { envelope_mass: m } => {
            assert_eq!(m, envelope_mass);
        }
        _ => panic!("Expected Hydrostatic envelope state"),
    }
}

// ===== Growth Physics Tests =====

#[test]
fn feeding_zone_scales_with_hill_radius() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let r_hill = body.hill_radius(Mass::from_solar_masses(1.0));
    let feeding_zone = body.feeding_zone_width(Mass::from_solar_masses(1.0));

    // Feeding zone should be 2.5 × Hill radius
    assert!(
        (feeding_zone.to_au() - 2.5 * r_hill.to_au()).abs() < 1e-10,
        "Feeding zone should be 2.5 × R_H"
    );
}

#[test]
fn isolation_mass_increases_with_surface_density() {
    use units::SurfaceDensity;

    let body = DiscreteBody::new(
        Mass::from_earth_masses(0.1),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma_low = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_high = SurfaceDensity::from_grams_per_cm2(100.0);

    let m_iso_low = body.isolation_mass(sigma_low, Mass::from_solar_masses(1.0));
    let m_iso_high = body.isolation_mass(sigma_high, Mass::from_solar_masses(1.0));

    assert!(
        m_iso_high > m_iso_low,
        "Higher surface density should give higher isolation mass"
    );

    // Should scale linearly with surface density
    let ratio = m_iso_high / m_iso_low;
    assert!(
        (ratio - 10.0).abs() < 0.1,
        "Isolation mass should scale linearly with Σ"
    );
}

#[test]
fn isolation_mass_scales_with_feeding_zone() {
    use units::SurfaceDensity;

    let small_body = DiscreteBody::new(
        Mass::from_earth_masses(0.1),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let large_body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma = SurfaceDensity::from_grams_per_cm2(100.0);

    let m_iso_small = small_body.isolation_mass(sigma, Mass::from_solar_masses(1.0));
    let m_iso_large = large_body.isolation_mass(sigma, Mass::from_solar_masses(1.0));

    // Larger body has wider feeding zone, thus higher isolation mass
    assert!(
        m_iso_large > m_iso_small,
        "Larger body should have higher isolation mass"
    );
}

#[test]
fn accretion_rate_positive() {
    use units::{SurfaceDensity, Velocity};

    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma_p = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_v = Velocity::from_meters_per_sec(100.0);

    let dm_dt = body.accretion_rate_from_particles(sigma_p, sigma_v, Mass::from_solar_masses(1.0));

    assert!(
        dm_dt.to_earth_masses_per_myr() > 0.0,
        "Accretion rate should be positive"
    );
}

#[test]
fn accretion_rate_increases_with_particle_density() {
    use units::{SurfaceDensity, Velocity};

    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma_low = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_high = SurfaceDensity::from_grams_per_cm2(100.0);
    let sigma_v = Velocity::from_meters_per_sec(100.0);

    let dm_dt_low =
        body.accretion_rate_from_particles(sigma_low, sigma_v, Mass::from_solar_masses(1.0));
    let dm_dt_high =
        body.accretion_rate_from_particles(sigma_high, sigma_v, Mass::from_solar_masses(1.0));

    assert!(
        dm_dt_high.to_earth_masses_per_myr() > dm_dt_low.to_earth_masses_per_myr(),
        "Higher particle density should give higher accretion rate"
    );

    // Should scale linearly with Σ_p
    let ratio = dm_dt_high.to_earth_masses_per_myr() / dm_dt_low.to_earth_masses_per_myr();
    assert!(
        (ratio - 10.0).abs() < 0.1,
        "Accretion rate should scale linearly with Σ_p"
    );
}

#[test]
fn accretion_rate_decreases_with_velocity_dispersion() {
    use units::{SurfaceDensity, Velocity};

    let body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let sigma_p = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_v_low = Velocity::from_meters_per_sec(10.0); // Low dispersion → strong focusing
    let sigma_v_high = Velocity::from_meters_per_sec(1000.0); // High dispersion → weak focusing

    let dm_dt_low =
        body.accretion_rate_from_particles(sigma_p, sigma_v_low, Mass::from_solar_masses(1.0));
    let dm_dt_high =
        body.accretion_rate_from_particles(sigma_p, sigma_v_high, Mass::from_solar_masses(1.0));

    assert!(
        dm_dt_low.to_earth_masses_per_myr() > dm_dt_high.to_earth_masses_per_myr(),
        "Lower velocity dispersion should give higher accretion rate (stronger focusing)"
    );
}

#[test]
fn accretion_rate_scales_with_radius_squared() {
    use units::{SurfaceDensity, Velocity};

    // Two bodies with same mass but different radii (due to composition)
    let small_body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(), // Dense, small radius
    );

    let fluffy_body = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        Mass::from_solar_masses(1.0),
        Composition::mini_neptune(), // Less dense, larger radius
    );

    let sigma_p = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_v = Velocity::from_meters_per_sec(100.0);

    let dm_dt_small =
        small_body.accretion_rate_from_particles(sigma_p, sigma_v, Mass::from_solar_masses(1.0));
    let dm_dt_fluffy =
        fluffy_body.accretion_rate_from_particles(sigma_p, sigma_v, Mass::from_solar_masses(1.0));

    // Larger radius → larger cross section → higher accretion rate
    assert!(
        dm_dt_fluffy.to_earth_masses_per_myr() > dm_dt_small.to_earth_masses_per_myr(),
        "Larger radius should give higher accretion rate"
    );
}

#[test]
fn runaway_growth_scaling_emerges() {
    use units::{SurfaceDensity, Velocity};

    // This is THE KEY TEST for runaway growth!
    // In the runaway regime (low velocity dispersion), accretion rate should
    // scale as dM/dt ∝ M^(4/3) due to gravitational focusing.
    //
    // dM/dt ∝ R² × F_g ∝ R² × (v_esc/σ)² ∝ R² × (M/R)²/σ² ∝ M²/σ²
    // But R ∝ M^(1/3) for constant density, so:
    // dM/dt ∝ M^(2/3) × (M^(2/3)/σ²) = M^(4/3)/σ²
    //
    // This is runaway growth - the bigger you are, the faster you grow!

    // Test with three bodies of different masses
    let masses = [0.1, 1.0, 10.0]; // Earth masses
    let stellar_mass = Mass::from_solar_masses(1.0);
    let sigma_p = SurfaceDensity::from_grams_per_cm2(10.0);
    let sigma_v = Velocity::from_meters_per_sec(50.0); // Low dispersion for runaway

    let mut accretion_rates = Vec::new();

    for &m_earth in &masses {
        let body = DiscreteBody::new(
            Mass::from_earth_masses(m_earth),
            Mass::zero(),
            Point2::new(1.0, 0.0),
            Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
            stellar_mass,
            Composition::earth_like(),
        );

        let dm_dt = body.accretion_rate_from_particles(sigma_p, sigma_v, stellar_mass);
        accretion_rates.push(dm_dt.to_earth_masses_per_myr());
    }

    // Check that dM/dt scales approximately as M^(4/3)
    // Compare ratios: (dM/dt)_2 / (dM/dt)_1 should ≈ (M_2 / M_1)^(4/3)

    let mass_ratio_1 = masses[1] / masses[0]; // 10.0
    let accretion_ratio_1 = accretion_rates[1] / accretion_rates[0];
    let expected_ratio_1 = mass_ratio_1.powf(4.0 / 3.0); // 10^(4/3) ≈ 21.5

    let mass_ratio_2 = masses[2] / masses[1]; // 10.0
    let accretion_ratio_2 = accretion_rates[2] / accretion_rates[1];
    let expected_ratio_2 = mass_ratio_2.powf(4.0 / 3.0);

    println!("Runaway Growth Test:");
    println!("  Mass ratio: {:.2}", mass_ratio_1);
    println!("  Expected accretion ratio: {:.2}", expected_ratio_1);
    println!("  Actual accretion ratio: {:.2}", accretion_ratio_1);
    println!();
    println!("  Mass ratio: {:.2}", mass_ratio_2);
    println!("  Expected accretion ratio: {:.2}", expected_ratio_2);
    println!("  Actual accretion ratio: {:.2}", accretion_ratio_2);

    // Allow 20% tolerance due to radius estimation and other factors
    assert!(
        (accretion_ratio_1 / expected_ratio_1 - 1.0).abs() < 0.2,
        "Accretion should scale as M^(4/3) in runaway regime. Expected {:.2}, got {:.2}",
        expected_ratio_1,
        accretion_ratio_1
    );

    assert!(
        (accretion_ratio_2 / expected_ratio_2 - 1.0).abs() < 0.2,
        "Accretion should scale as M^(4/3) in runaway regime. Expected {:.2}, got {:.2}",
        expected_ratio_2,
        accretion_ratio_2
    );

    // Also verify that larger bodies accrete faster (fundamental requirement)
    assert!(
        accretion_rates[1] > accretion_rates[0],
        "Runaway growth requires larger bodies to accrete faster"
    );
    assert!(
        accretion_rates[2] > accretion_rates[1],
        "Runaway growth requires larger bodies to accrete faster"
    );
}
