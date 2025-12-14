//! Simple orbital integration example
//!
//! Demonstrates the Leapfrog integrator with a single planet
//! orbiting a Sun-like star, showing energy conservation.
//!
//! Run with: cargo run --package nbody --example simple_orbit

use nalgebra::{Point2, Vector2};
use nbody::forces::{DirectGravity, ForceModel, G};
use nbody::integrator::{Integrator, Leapfrog};
use nbody::state::SystemState;
use stellar::generation::main_sequence_star;

fn main() {
    println!("N-body Integrator Test: Single Planet Orbit\n");
    println!("{}", "=".repeat(60));

    // Create a Sun-like star
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    println!("Star: M = {:.2} M☉", star.mass.to_solar_masses());

    // Create system state
    let mut system = SystemState::new(star);

    // Add Earth-like planet at 1 AU in circular orbit
    let r = 1.0; // AU
    let v_circular = (G * 1.0 / r).sqrt();
    let earth_mass = 3.0e-6; // M☉ (~ 1 Earth mass)
    let earth_radius = 4.26e-5; // AU (~ 1 Earth radius)

    system.add_body(
        earth_mass,
        earth_radius,
        Point2::new(r, 0.0),
        Vector2::new(0.0, v_circular),
    );

    println!("\nInitial conditions:");
    println!("  Planet mass: {:.2e} M☉", earth_mass);
    println!("  Orbital radius: {:.3} AU", r);
    println!("  Circular velocity: {:.3} AU/year", v_circular);

    // Set up integrator and force model
    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    // Initial diagnostics
    let initial_ke = system.kinetic_energy();
    let initial_pe = force.potential_energy(&system);
    let initial_energy = initial_ke + initial_pe;
    let initial_l = system.total_angular_momentum();

    println!("\nInitial energy:");
    println!("  Kinetic: {:.6e} M☉ AU² year⁻²", initial_ke);
    println!("  Potential: {:.6e} M☉ AU² year⁻²", initial_pe);
    println!("  Total: {:.6e} M☉ AU² year⁻²", initial_energy);
    println!("  Angular momentum: {:.6e} M☉ AU² year⁻¹", initial_l);

    // Integration parameters
    let orbital_period = 2.0 * std::f64::consts::PI; // 2π years for 1 AU
    let dt = orbital_period / 1000.0; // 1000 steps per orbit
    let n_orbits = 10;
    let n_steps = 1000 * n_orbits;

    println!("\nIntegration parameters:");
    println!("  Timestep: {:.6} years ({} steps/orbit)", dt, 1000);
    println!(
        "  Total time: {} orbits ({} years)",
        n_orbits,
        n_orbits as f64 * orbital_period
    );
    println!("  Total steps: {}", n_steps);

    println!("\nIntegrating...");

    // Integrate and track diagnostics
    let mut orbit_count = 0;
    let mut next_orbit_time = orbital_period;

    for _step in 0..n_steps {
        integrator.step(&mut system, dt, &force);

        // Print diagnostics every orbit
        if system.time >= next_orbit_time {
            orbit_count += 1;
            next_orbit_time += orbital_period;

            let body = &system.bodies[0];
            let r_current = body.position.coords.magnitude();
            let v_current = body.velocity.magnitude();

            let ke = system.kinetic_energy();
            let pe = force.potential_energy(&system);
            let energy = ke + pe;
            let l = system.total_angular_momentum();

            let energy_error = ((energy - initial_energy) / initial_energy).abs();
            let l_error = ((l - initial_l) / initial_l).abs();
            let r_error = ((r_current - r) / r).abs();

            println!(
                "Orbit {}: r={:.6} AU, v={:.4} AU/yr, ΔE={:.2e}, ΔL={:.2e}, Δr={:.2e}",
                orbit_count, r_current, v_current, energy_error, l_error, r_error
            );
        }
    }

    // Final diagnostics
    let final_ke = system.kinetic_energy();
    let final_pe = force.potential_energy(&system);
    let final_energy = final_ke + final_pe;
    let final_l = system.total_angular_momentum();

    println!("\n{}", "=".repeat(60));
    println!("Final diagnostics:");
    println!(
        "  Time: {:.2} years ({:.1} orbits)",
        system.time,
        system.time / orbital_period
    );

    let energy_error = ((final_energy - initial_energy) / initial_energy).abs();
    let l_error = ((final_l - initial_l) / initial_l).abs();

    println!("\nConservation:");
    println!(
        "  Energy error: {:.2e} ({:.4}%)",
        energy_error,
        energy_error * 100.0
    );
    println!(
        "  Angular momentum error: {:.2e} ({:.6}%)",
        l_error,
        l_error * 1e6
    );

    println!("\nFinal position:");
    let final_body = &system.bodies[0];
    println!(
        "  x = {:.6} AU, y = {:.6} AU",
        final_body.position.x, final_body.position.y
    );
    println!("  r = {:.6} AU", final_body.position.coords.magnitude());

    println!("\nFinal velocity:");
    println!(
        "  vx = {:.4} AU/yr, vy = {:.4} AU/yr",
        final_body.velocity.x, final_body.velocity.y
    );
    println!("  v = {:.4} AU/yr", final_body.velocity.magnitude());

    // Success criteria
    println!("\n{}", "=".repeat(60));
    if energy_error < 1e-4 {
        println!("✓ Energy conserved to within 0.01%");
    } else {
        println!("✗ Energy error too large: {:.2e}", energy_error);
    }

    if l_error < 1e-10 {
        println!("✓ Angular momentum conserved to machine precision");
    } else {
        println!("✗ Angular momentum error: {:.2e}", l_error);
    }

    let r_final = final_body.position.coords.magnitude();
    let r_error = ((r_final - r) / r).abs();
    if r_error < 0.01 {
        println!("✓ Orbit remains circular (Δr < 1%)");
    } else {
        println!("✗ Orbit drift: {:.2e}", r_error);
    }

    println!("\nTest complete!");
}
