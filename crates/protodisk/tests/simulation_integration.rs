//! Integration tests for the full simulation driver.
//!
//! These tests verify that all physics modules work together correctly.

use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use units::{Length, Mass, SurfaceDensity, Time};

use protodisk::bodies::DiscreteBody;
use protodisk::disk::{DiskMass, DiskModel, GasDisk, GridDisk};
use protodisk::particles::ParticleBin;
use protodisk::{SimulationState, run_simulation, solar_analog, step};

const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

#[test]
fn full_simulation_integration() {
    let star = solar_analog();
    let stellar_mass = star.mass;

    // Create gas disk
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

    // Create some particle bins at different radii
    let particle_bins = vec![
        ParticleBin::from_disk(&gas_disk, Length::from_au(1.0), Length::from_au(0.2)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(5.0), Length::from_au(0.5)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(10.0), Length::from_au(1.0)),
    ];

    // Create a proto-Jupiter at 5 AU
    let a = 5.0;
    let v_circ = (G * stellar_mass.to_solar_masses() / a).sqrt();
    let jupiter = DiscreteBody::new(
        Mass::from_earth_masses(10.0), // 10 M⊕ core
        Mass::zero(),                  // no envelope yet
        Point2::new(a, 0.0),
        Vector2::new(0.0, v_circ),
        stellar_mass,
        Composition::earth_like(),
    );

    // Create initial state
    let mut state = SimulationState::with_seed(grid_disk, particle_bins, vec![jupiter], 42);

    let initial_disk_mass = state.disk.total_mass();
    let initial_total_mass = state.total_mass();

    println!("\n=== Initial Conditions ===");
    println!("Time: {} years", state.time.to_years());
    println!("Disk mass: {} M☉", initial_disk_mass.to_solar_masses());
    println!("Total mass: {} M☉", initial_total_mass.to_solar_masses());
    println!("Particle bins: {}", state.particle_bins.len());
    println!("Discrete bodies: {}", state.discrete_bodies.len());
    println!(
        "Jupiter position: {:.2} AU, mass: {:.2} M⊕",
        state.discrete_bodies[0].semi_major_axis.to_au(),
        state.discrete_bodies[0].total_mass().to_earth_masses()
    );

    // Run for 100 timesteps
    println!("\n=== Running Simulation ===");
    for i in 0..100 {
        let dt = step(&mut state);

        if i % 10 == 0 {
            let jupiter = &state.discrete_bodies[0];
            println!(
                "Step {}: t={:.1} yr, dt={:.1} yr, Jupiter: {:.2} AU, {:.2} M⊕ (core: {:.2}, env: {:.2})",
                i,
                state.time.to_years(),
                dt.to_years(),
                jupiter.semi_major_axis.to_au(),
                jupiter.total_mass().to_earth_masses(),
                jupiter.core_mass.to_earth_masses(),
                jupiter.envelope_mass.to_earth_masses(),
            );
        }
    }

    println!("\n=== Final State ===");
    let final_disk_mass = state.disk.total_mass();
    let final_total_mass = state.total_mass();
    let jupiter = &state.discrete_bodies[0];

    println!("Time: {} years", state.time.to_years());
    println!("Disk mass: {} M☉", final_disk_mass.to_solar_masses());
    println!("Total mass: {} M☉", final_total_mass.to_solar_masses());
    println!(
        "Jupiter: {:.2} AU, {:.2} M⊕ (core: {:.2}, env: {:.2})",
        jupiter.semi_major_axis.to_au(),
        jupiter.total_mass().to_earth_masses(),
        jupiter.core_mass.to_earth_masses(),
        jupiter.envelope_mass.to_earth_masses(),
    );

    // Verify physics worked
    println!("\n=== Physics Validation ===");

    // 1. Disk should have lost mass (photoevaporation)
    assert!(
        final_disk_mass < initial_disk_mass,
        "Disk should lose mass to photoevaporation"
    );
    println!("✓ Disk evolution: Mass decreased from photoevaporation");

    // 2. Jupiter should have migrated inward
    assert!(
        jupiter.semi_major_axis.to_au() < a,
        "Jupiter should migrate inward"
    );
    println!(
        "✓ Migration: Jupiter migrated inward {:.2} AU",
        a - jupiter.semi_major_axis.to_au()
    );

    // 3. Jupiter should have captured an envelope
    assert!(
        jupiter.envelope_mass.to_earth_masses() > 0.0,
        "Jupiter should have gas envelope"
    );
    println!(
        "✓ Envelope evolution: Captured {:.2} M⊕ of gas",
        jupiter.envelope_mass.to_earth_masses()
    );

    // 4. Mass should be approximately conserved (except photoevaporation losses)
    let mass_lost = initial_total_mass - final_total_mass;
    let fractional_loss = mass_lost.to_solar_masses() / initial_total_mass.to_solar_masses();
    assert!(
        fractional_loss < 0.1,
        "Should not lose more than 10% of mass"
    );
    println!(
        "✓ Mass conservation: {:.4}% lost (expected from photoevaporation)",
        fractional_loss * 100.0
    );

    // 5. Time should have advanced
    assert!(state.time.to_years() > 0.0);
    println!(
        "✓ Time integration: Advanced {:.1} years",
        state.time.to_years()
    );

    println!("\n=== All Physics Tests Passed ===");
}

#[test]
fn planetesimal_formation_from_unstable_layer() {
    let star = solar_analog();
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

    // Create a dense, settled particle bin (gravitationally unstable)
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);
    let mut bin = ParticleBin::from_disk(&gas_disk, r, width);

    // Make it very dense and settled (unstable conditions)
    // Use same parameters as particle_bin_test::unstable_bin_forms_planetesimals
    bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(100.0));
    bin.set_scale_height(grid_disk.scale_height(r) * 0.5);
    bin.set_velocity_dispersion(grid_disk.sound_speed(r) * 0.001);

    let initial_mass = bin.total_mass();

    println!("\n=== Planetesimal Formation Test ===");
    println!(
        "Initial bin mass: {:.2e} M⊕",
        initial_mass.to_earth_masses()
    );
    println!("Toomre Q: {:.3}", bin.toomre_q(&gas_disk));
    println!("Richardson Ri: {:.3}", bin.richardson_number(&gas_disk));
    println!(
        "Gravitationally unstable: {}",
        bin.is_gravitationally_unstable(&gas_disk)
    );

    let mut state = SimulationState::with_seed(grid_disk, vec![bin], vec![], 42);

    // Run for a few steps
    let mut total_formed_mass = Mass::zero();
    for i in 0..10 {
        step(&mut state);

        if i % 2 == 0 {
            let bin_mass = if !state.particle_bins.is_empty() {
                state.particle_bins[0].total_mass()
            } else {
                Mass::zero()
            };

            let formed = initial_mass - bin_mass;
            total_formed_mass = formed;

            println!(
                "Step {}: Bin mass: {:.2e} M⊕, Formed: {:.2e} M⊕",
                i,
                bin_mass.to_earth_masses(),
                formed.to_earth_masses()
            );
        }
    }

    println!("\n=== Results ===");
    println!(
        "Total planetesimal mass formed: {:.2e} M⊕",
        total_formed_mass.to_earth_masses()
    );

    // Should have formed some planetesimals
    assert!(
        total_formed_mass.to_earth_masses() > 0.0,
        "Should form planetesimals from unstable layer"
    );
    println!("✓ Planetesimal formation working");
}

#[test]
fn multiple_planets_evolve_independently() {
    let star = solar_analog();
    let stellar_mass = star.mass;
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

    // Create three protoplanets at different locations
    let planets = vec![
        // Inner super-Earth at 1 AU
        {
            let a = 1.0;
            let v = (G * stellar_mass.to_solar_masses() / a).sqrt();
            DiscreteBody::new(
                Mass::from_earth_masses(5.0),
                Mass::zero(),
                Point2::new(a, 0.0),
                Vector2::new(0.0, v),
                stellar_mass,
                Composition::earth_like(),
            )
        },
        // Jupiter-like at 5 AU
        {
            let a = 5.0;
            let v = (G * stellar_mass.to_solar_masses() / a).sqrt();
            DiscreteBody::new(
                Mass::from_earth_masses(10.0),
                Mass::zero(),
                Point2::new(a, 0.0),
                Vector2::new(0.0, v),
                stellar_mass,
                Composition::earth_like(),
            )
        },
        // Outer ice giant at 10 AU
        {
            let a = 10.0;
            let v = (G * stellar_mass.to_solar_masses() / a).sqrt();
            DiscreteBody::new(
                Mass::from_earth_masses(8.0),
                Mass::zero(),
                Point2::new(a, 0.0),
                Vector2::new(0.0, v),
                stellar_mass,
                Composition::earth_like(),
            )
        },
    ];

    let mut state = SimulationState::with_seed(grid_disk, vec![], planets, 42);

    println!("\n=== Multiple Planets Test ===");
    println!("Initial configuration:");
    for (i, planet) in state.discrete_bodies.iter().enumerate() {
        println!(
            "  Planet {}: {:.1} AU, {:.1} M⊕",
            i,
            planet.semi_major_axis.to_au(),
            planet.total_mass().to_earth_masses()
        );
    }

    // Run simulation
    for _ in 0..50 {
        step(&mut state);
    }

    println!("\nFinal configuration:");
    for (i, planet) in state.discrete_bodies.iter().enumerate() {
        println!(
            "  Planet {}: {:.1} AU, {:.1} M⊕ (env: {:.1})",
            i,
            planet.semi_major_axis.to_au(),
            planet.total_mass().to_earth_masses(),
            planet.envelope_mass.to_earth_masses()
        );
    }

    // Verify all planets evolved
    assert_eq!(state.discrete_bodies.len(), 3);

    // Inner planet should have migrated inward
    assert!(state.discrete_bodies[0].semi_major_axis.to_au() < 1.0);

    // Middle planet should have captured envelope
    assert!(state.discrete_bodies[1].envelope_mass.to_earth_masses() > 0.0);

    // Outer planet should have migrated
    assert!(state.discrete_bodies[2].semi_major_axis.to_au() < 10.0);

    println!("✓ All planets evolved independently");
}

#[test]
fn run_until_completion() {
    let star = solar_analog();
    let stellar_mass = star.mass;
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

    // Create a super-Earth
    let a = 3.0;
    let v = (G * stellar_mass.to_solar_masses() / a).sqrt();
    let planet = DiscreteBody::new(
        Mass::from_earth_masses(8.0),
        Mass::zero(),
        Point2::new(a, 0.0),
        Vector2::new(0.0, v),
        stellar_mass,
        Composition::earth_like(),
    );

    let state = SimulationState::with_seed(grid_disk, vec![], vec![planet], 42);

    println!("\n=== Run to Completion Test ===");
    println!("Running simulation for up to 1 Myr...");

    let max_time = Time::from_years(1_000_000.0);
    let final_state = run_simulation(state, max_time);

    println!("Simulation completed!");
    println!("Final time: {:.2e} years", final_state.time.to_years());
    println!(
        "Disk mass: {:.2e} M☉",
        final_state.disk.total_mass().to_solar_masses()
    );
    println!("Bodies: {}", final_state.discrete_bodies.len());

    if let Some(planet) = final_state.discrete_bodies.first() {
        println!(
            "Planet: {:.2} AU, {:.2} M⊕",
            planet.semi_major_axis.to_au(),
            planet.total_mass().to_earth_masses()
        );
    }

    // Simulation should have finished
    assert!(final_state.is_finished(max_time));
    println!("✓ Simulation ran to completion");
}
