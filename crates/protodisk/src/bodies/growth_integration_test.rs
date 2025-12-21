//! Integration tests for planet growth from particle accretion.
//!
//! These tests validate that the growth physics works end-to-end:
//! - Particles exist in bins with surface density
//! - Discrete bodies accrete from those bins
//! - Mass flows from particles → bodies over time
//! - Particle bins deplete as bodies grow

use super::*;
use crate::disk::GasDisk;
use crate::particles::ParticleBin;
use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use units::{Length, Mass, SurfaceDensity};

const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

#[test]
fn planetesimal_grows_by_accreting_particles() {
    // This is the KEY integration test!
    // We create a particle bin and a small body, then watch it grow.

    println!("\n=== Planetesimal Growth Integration Test ===\n");

    // 1. Set up the environment at 1 AU
    let disk = GasDisk::mmsn();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);
    let stellar_mass = Mass::from_solar_masses(1.0);

    // 2. Create a particle bin with substantial surface density
    let mut particle_bin = ParticleBin::from_disk(&disk, r, width);

    // Boost particle surface density to speed up growth for testing
    let initial_particle_sigma = SurfaceDensity::from_grams_per_cm2(100.0);
    particle_bin.set_surface_density(initial_particle_sigma);

    let initial_particle_mass = particle_bin.total_mass();
    println!(
        "Initial particle mass in bin: {:.3} M⊕",
        initial_particle_mass.to_earth_masses()
    );

    // 3. Create a small discrete body (Moon-mass planetesimal)
    let initial_body_mass = Mass::from_earth_masses(0.01); // ~1 Moon mass
    let mut body = DiscreteBody::new(
        initial_body_mass,
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        stellar_mass,
        Composition::earth_like(),
    );

    println!(
        "Initial body mass: {:.3} M⊕",
        initial_body_mass.to_earth_masses()
    );
    println!(
        "Body radius: {:.4} R⊕",
        body.physical_radius.to_earth_radii()
    );
    println!(
        "Hill radius: {:.4} AU",
        body.hill_radius(stellar_mass).to_au()
    );

    // 4. Run growth simulation for 100,000 years (100 steps of 1000 years)
    let timestep = 1000.0; // years
    let num_steps = 100;

    println!(
        "\nSimulating growth over {} timesteps of {} years...\n",
        num_steps, timestep
    );

    for step in 0..num_steps {
        // Calculate current accretion rate
        let particle_sigma = particle_bin.surface_density();
        let particle_sigma_v = particle_bin.velocity_dispersion();

        let dm_dt =
            body.accretion_rate_from_particles(particle_sigma, particle_sigma_v, stellar_mass);

        // Apply growth: ΔM = (dM/dt) × Δt
        let mass_wanted = dm_dt.to_solar_masses_per_year() * timestep;
        let mass_wanted_mass = Mass::from_solar_masses(mass_wanted);

        // Limit accretion to available mass in particle bin
        let mass_available = particle_bin.total_mass();
        let mass_accreted_mass = if mass_wanted_mass > mass_available {
            mass_available * 0.9 // Take at most 90% to avoid complete depletion
        } else {
            mass_wanted_mass
        };

        // Update body mass
        body.core_mass = body.core_mass + mass_accreted_mass;

        // Update body radius (mass changed)
        body.physical_radius = DiscreteBody::estimate_radius(body.total_mass(), &body.composition);

        // Remove mass from particle bin
        let sigma_to_remove = mass_accreted_mass.to_grams() / particle_bin.area();
        let new_particle_sigma = SurfaceDensity::from_grams_per_cm2(
            (particle_sigma.to_grams_per_cm2() - sigma_to_remove).max(0.0),
        );
        particle_bin.set_surface_density(new_particle_sigma);

        // Log every 20 steps
        if step % 20 == 0 {
            let time_kyr = (step as f64 * timestep) / 1000.0;
            println!(
                "t = {:.0} kyr: M = {:.4} M⊕, R = {:.4} R⊕, dM/dt = {:.2} M⊕/Myr, Σ_p = {:.2} g/cm²",
                time_kyr,
                body.total_mass().to_earth_masses(),
                body.physical_radius.to_earth_radii(),
                dm_dt.to_earth_masses_per_myr(),
                particle_sigma.to_grams_per_cm2(),
            );
        }
    }

    let final_body_mass = body.total_mass();
    let final_particle_mass = particle_bin.total_mass();

    println!("\n=== Final State ===");
    println!(
        "Body mass: {:.4} M⊕ (grew {:.2}×)",
        final_body_mass.to_earth_masses(),
        final_body_mass / initial_body_mass
    );
    println!(
        "Particle mass: {:.3} M⊕ (depleted to {:.1}%)",
        final_particle_mass.to_earth_masses(),
        100.0 * final_particle_mass / initial_particle_mass
    );

    // Verify mass conservation
    let total_initial = initial_body_mass + initial_particle_mass;
    let total_final = final_body_mass + final_particle_mass;

    println!("\nMass conservation check:");
    println!("  Initial total: {:.4} M⊕", total_initial.to_earth_masses());
    println!("  Final total:   {:.4} M⊕", total_final.to_earth_masses());
    println!(
        "  Difference:    {:.2}%",
        100.0 * (total_final / total_initial - 1.0)
    );

    // Assertions
    assert!(
        final_body_mass > initial_body_mass * 2.0,
        "Body should have grown significantly (at least doubled)"
    );

    assert!(
        final_particle_mass < initial_particle_mass * 0.95,
        "Particle bin should have depleted"
    );

    // Mass should be conserved to within numerical precision
    assert!(
        (total_final / total_initial - 1.0).abs() < 1e-6,
        "Mass should be conserved"
    );

    println!("\n✓ All checks passed! Planet growth works!\n");
}

#[test]
fn runaway_growth_accelerates_over_time() {
    // Verify that accretion rate INCREASES as the body grows (runaway!)

    println!("\n=== Runaway Growth Acceleration Test ===\n");

    let disk = GasDisk::mmsn();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);
    let stellar_mass = Mass::from_solar_masses(1.0);

    let mut particle_bin = ParticleBin::from_disk(&disk, r, width);
    particle_bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(500.0)); // Much larger reservoir

    // Start with a small body
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(0.01), // Smaller starting mass
        Mass::zero(),
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
        stellar_mass,
        Composition::earth_like(),
    );

    let timestep = 500.0; // years (shorter timesteps)
    let mut accretion_rates = Vec::new();
    let mut masses = Vec::new();

    // Track growth for 30 timesteps (before depletion)
    for _ in 0..30 {
        let dm_dt = body.accretion_rate_from_particles(
            particle_bin.surface_density(),
            particle_bin.velocity_dispersion(),
            stellar_mass,
        );

        accretion_rates.push(dm_dt.to_earth_masses_per_myr());
        masses.push(body.total_mass().to_earth_masses());

        // Apply growth (limit to available mass)
        let mass_wanted = Mass::from_solar_masses(dm_dt.to_solar_masses_per_year() * timestep);
        let mass_available = particle_bin.total_mass();
        let mass_accreted = if mass_wanted > mass_available {
            mass_available * 0.9
        } else {
            mass_wanted
        };

        body.core_mass = body.core_mass + mass_accreted;
        body.physical_radius = DiscreteBody::estimate_radius(body.total_mass(), &body.composition);

        // Deplete particles
        let sigma_to_remove = mass_accreted.to_grams() / particle_bin.area();
        let new_sigma = SurfaceDensity::from_grams_per_cm2(
            (particle_bin.surface_density().to_grams_per_cm2() - sigma_to_remove).max(0.0),
        );
        particle_bin.set_surface_density(new_sigma);
    }

    println!("Initial accretion rate: {:.4} M⊕/Myr", accretion_rates[0]);
    println!(
        "Peak accretion rate:    {:.4} M⊕/Myr",
        accretion_rates.iter().cloned().fold(0.0f64, f64::max)
    );
    println!(
        "Final accretion rate:   {:.4} M⊕/Myr",
        accretion_rates[accretion_rates.len() - 1]
    );

    println!("\nInitial mass: {:.4} M⊕", masses[0]);
    println!("Final mass:   {:.4} M⊕", masses[masses.len() - 1]);
    println!(
        "Growth factor: {:.2}×",
        masses[masses.len() - 1] / masses[0]
    );

    // In runaway regime, there should be a period where accretion accelerates
    // Check that the peak accretion rate is higher than initial (before depletion kicks in)
    let peak_rate = accretion_rates.iter().cloned().fold(0.0f64, f64::max);
    assert!(
        peak_rate > accretion_rates[0] * 1.2,
        "Peak accretion rate should be higher than initial (runaway growth)"
    );

    // Find where accretion peaked and verify it was due to mass increase
    let peak_idx = accretion_rates
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    assert!(
        masses[peak_idx] > masses[0],
        "Peak accretion should occur at higher mass (runaway regime)"
    );

    println!(
        "\n✓ Runaway growth confirmed! Accretion accelerated from {:.2} to {:.2} M⊕/Myr as mass increased.\n",
        accretion_rates[0], peak_rate
    );
}

#[test]
fn multiple_bodies_deplete_shared_particle_bin() {
    // Test that multiple protoplanets can accrete from the same particle reservoir

    println!("\n=== Multi-Body Competition Test ===\n");

    let disk = GasDisk::mmsn();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.2); // Wider bin
    let stellar_mass = Mass::from_solar_masses(1.0);

    let mut particle_bin = ParticleBin::from_disk(&disk, r, width);
    particle_bin.set_surface_density(SurfaceDensity::from_grams_per_cm2(200.0));

    let initial_particle_mass = particle_bin.total_mass();

    // Create 3 competing embryos
    let mut bodies = vec![
        DiscreteBody::new(
            Mass::from_earth_masses(0.1),
            Mass::zero(),
            Point2::new(0.95, 0.0),
            Vector2::new(0.0, (G * 1.0 / 0.95).sqrt()),
            stellar_mass,
            Composition::earth_like(),
        ),
        DiscreteBody::new(
            Mass::from_earth_masses(0.1),
            Mass::zero(),
            Point2::new(1.0, 0.0),
            Vector2::new(0.0, (G * 1.0 / 1.0).sqrt()),
            stellar_mass,
            Composition::earth_like(),
        ),
        DiscreteBody::new(
            Mass::from_earth_masses(0.1),
            Mass::zero(),
            Point2::new(1.05, 0.0),
            Vector2::new(0.0, (G * 1.0 / 1.05).sqrt()),
            stellar_mass,
            Composition::earth_like(),
        ),
    ];

    println!("Created 3 embryos of 0.1 M⊕ each");
    println!(
        "Particle bin: {:.2} M⊕\n",
        initial_particle_mass.to_earth_masses()
    );

    let timestep = 500.0; // years
    let num_steps = 50;

    for step in 0..num_steps {
        let mut total_accreted = Mass::zero();

        // Each body accretes
        let mass_available = particle_bin.total_mass();

        for body in &mut bodies {
            let dm_dt = body.accretion_rate_from_particles(
                particle_bin.surface_density(),
                particle_bin.velocity_dispersion(),
                stellar_mass,
            );

            let mass_wanted = Mass::from_solar_masses(dm_dt.to_solar_masses_per_year() * timestep);

            // Limit each body to fair share of remaining mass
            let mass_accreted = if total_accreted + mass_wanted > mass_available {
                (mass_available - total_accreted) * 0.3 // Take small fraction of what's left
            } else {
                mass_wanted
            };

            body.core_mass = body.core_mass + mass_accreted;
            body.physical_radius =
                DiscreteBody::estimate_radius(body.total_mass(), &body.composition);

            total_accreted = total_accreted + mass_accreted;
        }

        // Deplete particle bin by total accreted
        let sigma_to_remove = total_accreted.to_grams() / particle_bin.area();
        let new_sigma = SurfaceDensity::from_grams_per_cm2(
            (particle_bin.surface_density().to_grams_per_cm2() - sigma_to_remove).max(0.0),
        );
        particle_bin.set_surface_density(new_sigma);

        if step % 10 == 0 {
            let avg_mass = bodies
                .iter()
                .map(|b| b.total_mass().to_earth_masses())
                .sum::<f64>()
                / bodies.len() as f64;

            println!(
                "t = {:.0} kyr: Avg mass = {:.4} M⊕, Particles = {:.2} M⊕",
                (step as f64 * timestep) / 1000.0,
                avg_mass,
                particle_bin.total_mass().to_earth_masses()
            );
        }
    }

    let final_particle_mass = particle_bin.total_mass();
    let total_body_mass: Mass = bodies
        .iter()
        .map(|b| b.total_mass())
        .fold(Mass::zero(), |acc, m| acc + m);

    println!("\n=== Final State ===");
    for (i, body) in bodies.iter().enumerate() {
        println!(
            "Body {}: {:.4} M⊕",
            i + 1,
            body.total_mass().to_earth_masses()
        );
    }
    println!("Particles: {:.2} M⊕", final_particle_mass.to_earth_masses());
    println!(
        "Total in bodies: {:.4} M⊕",
        total_body_mass.to_earth_masses()
    );

    // Verify all bodies grew
    for body in &bodies {
        assert!(
            body.total_mass() > Mass::from_earth_masses(0.15),
            "Each body should have grown significantly"
        );
    }

    // Verify particles depleted
    assert!(
        final_particle_mass < initial_particle_mass * 0.9,
        "Particle bin should be significantly depleted"
    );

    println!("\n✓ Multi-body competition works correctly!\n");
}
