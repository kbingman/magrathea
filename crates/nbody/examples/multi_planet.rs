//! Multi-planet system integration
//!
//! Demonstrates N-body integration with multiple interacting planets,
//! similar to an inner solar system.
//!
//! Run with: cargo run --package nbody --example multi_planet

use nalgebra::{Point2, Vector2};
use nbody::forces::{DirectGravity, ForceModel, G, TreeGravity};
use nbody::integrator::{Integrator, Leapfrog};
use nbody::state::SystemState;
use stellar::generation::main_sequence_star;
use units::{Length, Mass};

fn main() {
    println!("N-body Integrator Test: Multi-Planet System\n");
    println!("{}", "=".repeat(60));

    // Create a Sun-like star
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add inner solar system analogs
    // Mercury: 0.39 AU, 0.055 M⊕
    let mercury_mass = Mass::from_earth_masses(0.055).to_solar_masses();
    let mercury_r = 0.39;
    let mercury_v = (G * 1.0 / mercury_r).sqrt();
    system.add_body(
        mercury_mass,
        Length::from_earth_radii(0.38).to_au(),
        Point2::new(mercury_r, 0.0),
        Vector2::new(0.0, mercury_v),
    );

    // Venus: 0.72 AU, 0.815 M⊕
    let venus_mass = Mass::from_earth_masses(0.815).to_solar_masses();
    let venus_r = 0.72;
    let venus_v = (G * 1.0 / venus_r).sqrt();
    system.add_body(
        venus_mass,
        Length::from_earth_radii(0.95).to_au(),
        Point2::new(venus_r, 0.0),
        Vector2::new(0.0, venus_v),
    );

    // Earth: 1.0 AU, 1.0 M⊕
    let earth_mass = Mass::from_earth_masses(1.0).to_solar_masses();
    let earth_r = 1.0;
    let earth_v = (G * 1.0 / earth_r).sqrt();
    system.add_body(
        earth_mass,
        Length::from_earth_radii(1.0).to_au(),
        Point2::new(earth_r, 0.0),
        Vector2::new(0.0, earth_v),
    );

    // Mars: 1.52 AU, 0.107 M⊕
    let mars_mass = Mass::from_earth_masses(0.107).to_solar_masses();
    let mars_r = 1.52;
    let mars_v = (G * 1.0 / mars_r).sqrt();
    system.add_body(
        mars_mass,
        Length::from_earth_radii(0.53).to_au(),
        Point2::new(mars_r, 0.0),
        Vector2::new(0.0, mars_v),
    );

    println!("\nInitial system:");
    println!("  Number of planets: {}", system.body_count());
    println!("  Total planet mass: {:.3e} M☉", system.total_planet_mass());
    println!(
        "  Total planet mass: {:.3} M⊕",
        Mass::from_solar_masses(system.total_planet_mass()).to_earth_masses()
    );

    // Test both force models
    let force_direct = DirectGravity::new();
    let force_tree = TreeGravity::new();

    let initial_ke = system.kinetic_energy();
    let initial_pe_direct = force_direct.potential_energy(&system);
    let initial_pe_tree = force_tree.potential_energy(&system);
    let initial_l = system.total_angular_momentum();

    println!("\nInitial energy (direct gravity):");
    println!("  Kinetic: {:.6e} M☉ AU² year⁻²", initial_ke);
    println!("  Potential: {:.6e} M☉ AU² year⁻²", initial_pe_direct);
    println!(
        "  Total: {:.6e} M☉ AU² year⁻²",
        initial_ke + initial_pe_direct
    );

    println!("\nInitial energy (tree gravity):");
    println!("  Potential: {:.6e} M☉ AU² year⁻²", initial_pe_tree);
    println!(
        "  Total: {:.6e} M☉ AU² year⁻²",
        initial_ke + initial_pe_tree
    );
    println!(
        "  Difference: {:.2e}",
        (initial_pe_tree - initial_pe_direct).abs()
    );

    println!("\nAngular momentum:");
    println!("  Total: {:.6e} M☉ AU² year⁻¹", initial_l);

    // Integration parameters
    let dt = 0.01; // years
    let n_years = 100;
    let n_steps = (n_years as f64 / dt) as usize;

    println!("\nIntegration parameters:");
    println!("  Timestep: {:.3} years", dt);
    println!("  Total time: {} years", n_years);
    println!("  Total steps: {}", n_steps);

    // Clone system for tree comparison
    let mut system_tree = system.clone();

    println!("\nIntegrating with DirectGravity (O(N²))...");
    let integrator = Leapfrog::new();

    let mut print_interval = 10.0;
    let mut next_print = print_interval;

    for _step in 0..n_steps {
        integrator.step(&mut system, dt, &force_direct);

        if system.time >= next_print {
            let ke = system.kinetic_energy();
            let pe = force_direct.potential_energy(&system);
            let energy = ke + pe;
            let l = system.total_angular_momentum();

            let energy_error = ((energy - (initial_ke + initial_pe_direct))
                / (initial_ke + initial_pe_direct))
                .abs();
            let l_error = ((l - initial_l) / initial_l).abs();

            println!(
                "  t={:5.0} yr: ΔE={:.2e}, ΔL={:.2e}",
                system.time, energy_error, l_error
            );

            next_print += print_interval;
        }
    }

    println!("\nIntegrating with TreeGravity (O(N log N))...");

    print_interval = 10.0;
    next_print = print_interval;

    for _step in 0..n_steps {
        integrator.step(&mut system_tree, dt, &force_tree);

        if system_tree.time >= next_print {
            let ke = system_tree.kinetic_energy();
            let pe = force_tree.potential_energy(&system_tree);
            let energy = ke + pe;
            let l = system_tree.total_angular_momentum();

            let energy_error =
                ((energy - (initial_ke + initial_pe_tree)) / (initial_ke + initial_pe_tree)).abs();
            let l_error = ((l - initial_l) / initial_l).abs();

            println!(
                "  t={:5.0} yr: ΔE={:.2e}, ΔL={:.2e}",
                system_tree.time, energy_error, l_error
            );

            next_print += print_interval;
        }
    }

    // Final comparison
    println!("\n{}", "=".repeat(60));
    println!("Final results after {} years:", n_years);

    println!("\nDirectGravity:");
    let final_ke_direct = system.kinetic_energy();
    let final_pe_direct = force_direct.potential_energy(&system);
    let final_energy_direct = final_ke_direct + final_pe_direct;
    let final_l_direct = system.total_angular_momentum();

    let energy_error_direct = ((final_energy_direct - (initial_ke + initial_pe_direct))
        / (initial_ke + initial_pe_direct))
        .abs();
    let l_error_direct = ((final_l_direct - initial_l) / initial_l).abs();

    println!(
        "  Energy error: {:.2e} ({:.4}%)",
        energy_error_direct,
        energy_error_direct * 100.0
    );
    println!("  Angular momentum error: {:.2e}", l_error_direct);

    println!("\nTreeGravity:");
    let final_ke_tree = system_tree.kinetic_energy();
    let final_pe_tree = force_tree.potential_energy(&system_tree);
    let final_energy_tree = final_ke_tree + final_pe_tree;
    let final_l_tree = system_tree.total_angular_momentum();

    let energy_error_tree = ((final_energy_tree - (initial_ke + initial_pe_tree))
        / (initial_ke + initial_pe_tree))
        .abs();
    let l_error_tree = ((final_l_tree - initial_l) / initial_l).abs();

    println!(
        "  Energy error: {:.2e} ({:.4}%)",
        energy_error_tree,
        energy_error_tree * 100.0
    );
    println!("  Angular momentum error: {:.2e}", l_error_tree);

    // Compare final positions
    println!("\nFinal planet positions (DirectGravity):");
    let names = ["Mercury", "Venus", "Earth", "Mars"];
    for (i, name) in names.iter().enumerate() {
        let body = &system.bodies[i];
        let r = body.position.coords.magnitude();
        let v = body.velocity.magnitude();
        println!("  {}: r={:.4} AU, v={:.3} AU/yr", name, r, v);
    }

    println!("\nPosition differences (TreeGravity vs DirectGravity):");
    for (i, name) in names.iter().enumerate() {
        let pos_diff = (system_tree.bodies[i].position - system.bodies[i].position).magnitude();
        let vel_diff = (system_tree.bodies[i].velocity - system.bodies[i].velocity).magnitude();
        println!(
            "  {}: Δr={:.2e} AU, Δv={:.2e} AU/yr",
            name, pos_diff, vel_diff
        );
    }

    println!("\n{}", "=".repeat(60));
    println!("Test complete!");
}
