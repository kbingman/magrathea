//! Collision detection and resolution example
//!
//! Demonstrates collision handling with multiple bodies on intersecting orbits.
//! Shows how bodies merge when they come within their collision radius.
//!
//! Run with: cargo run --package nbody --example collision_demo

use nalgebra::{Point2, Vector2};
use nbody::collisions::{
    CollisionCriteria, CollisionDetector, DirectDetector, detect_ejections, detect_star_collisions,
    remove_ejections, remove_star_collisions, resolve_collisions,
};
use nbody::forces::{DirectGravity, G};
use nbody::integrator::{Integrator, Leapfrog};
use nbody::state::SystemState;
use stellar::generation::main_sequence_star;
use units::Mass;

fn main() {
    println!("N-body Collision Demo: Giant Impact Simulation\n");
    println!("{}", "=".repeat(60));

    // Create a Sun-like star
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add several planetary embryos on slightly different orbits
    // This simulates the late stage of planet formation with crossing orbits
    let earth_mass = Mass::from_earth_masses(1.0).to_solar_masses();

    println!("\nInitial system:");
    println!("  Creating 8 planetary embryos on crossing orbits...\n");

    for i in 0..8 {
        let a = 0.9 + (i as f64) * 0.025; // Semi-major axes from 0.9 to 1.075 AU
        let e = 0.02 + (i as f64) * 0.005; // Slight eccentricity variation
        let angle = (i as f64) * std::f64::consts::PI / 4.0; // Spread around orbit

        // Start at different positions
        let x = a * angle.cos();
        let y = a * angle.sin();

        // Circular orbital velocity
        let v_mag = (G * 1.0 / a).sqrt();
        let vx = -v_mag * angle.sin();
        let vy = v_mag * angle.cos();

        // Perturb slightly to create crossing orbits
        let vx_pert = vx * (1.0 + e);
        let vy_pert = vy * (1.0 + e);

        system.add_body(
            earth_mass * 0.5, // Half Earth mass each
            0.005,            // 0.005 AU radius
            Point2::new(x, y),
            Vector2::new(vx_pert, vy_pert),
        );

        println!(
            "  Embryo {}: a={:.3} AU, e={:.3}, pos=({:.3}, {:.3}) AU",
            i + 1,
            a,
            e,
            x,
            y
        );
    }

    println!("\nInitial planet count: {}", system.body_count());
    println!(
        "Initial total mass: {:.3} M⊕",
        Mass::from_solar_masses(system.total_planet_mass()).to_earth_masses()
    );

    // Integration setup
    let integrator = Leapfrog::new();
    let force = DirectGravity::new();
    let detector = DirectDetector;

    let criteria = CollisionCriteria {
        hill_fraction: 0.5, // Half Hill radius
        physical_collision: true,
    };

    let dt = 0.001; // Small timestep for accuracy
    let total_years = 10.0;
    let steps_per_check = 100; // Check for collisions every 100 steps
    let total_steps = (total_years / dt) as usize;

    println!("\nSimulation parameters:");
    println!("  Timestep: {} years", dt);
    println!("  Total time: {} years", total_years);
    println!("  Collision check interval: {} steps", steps_per_check);

    println!("\n{}", "=".repeat(60));
    println!("Starting simulation...\n");

    let mut collision_count = 0;
    let mut ejection_count = 0;
    let mut star_collision_count = 0;

    for step in 0..total_steps {
        integrator.step(&mut system, dt, &force);

        // Check for collisions periodically
        if step % steps_per_check == 0 {
            // Detect and resolve collisions
            let collisions = detector.detect(&system, &criteria);
            if !collisions.is_empty() {
                println!(
                    "t={:6.3} yr: {} collision(s) detected!",
                    system.time,
                    collisions.len()
                );

                for event in &collisions {
                    println!(
                        "    Bodies {:?} and {:?}: separation={:.4} AU, threshold={:.4} AU",
                        event.body_a, event.body_b, event.separation, event.collision_radius
                    );
                }

                collision_count += collisions.len();
                resolve_collisions(&mut system, collisions);

                println!("    Resolved! Bodies merged.");
                println!(
                    "    New planet count: {}, Total mass: {:.3} M⊕\n",
                    system.body_count(),
                    Mass::from_solar_masses(system.total_planet_mass()).to_earth_masses()
                );
            }

            // Check for ejections
            let ejections = detect_ejections(&system, 5.0); // Beyond 5 AU
            if !ejections.is_empty() {
                println!(
                    "t={:6.3} yr: {} body(ies) ejected!",
                    system.time,
                    ejections.len()
                );
                ejection_count += ejections.len();
                remove_ejections(&mut system, ejections);
            }

            // Check for star collisions
            let star_collisions = detect_star_collisions(&system);
            if !star_collisions.is_empty() {
                println!(
                    "t={:6.3} yr: {} body(ies) fell into star!",
                    system.time,
                    star_collisions.len()
                );
                star_collision_count += star_collisions.len();
                remove_star_collisions(&mut system, star_collisions);
            }
        }
    }

    // Final results
    println!("{}", "=".repeat(60));
    println!("Simulation complete!\n");

    println!("Final statistics:");
    println!("  Final time: {:.2} years", system.time);
    println!("  Final planet count: {}", system.body_count());
    println!(
        "  Final total mass: {:.3} M⊕",
        Mass::from_solar_masses(system.total_planet_mass()).to_earth_masses()
    );
    println!("  Total collisions: {}", collision_count);
    println!("  Total ejections: {}", ejection_count);
    println!("  Total star collisions: {}", star_collision_count);

    if !system.bodies.is_empty() {
        println!("\nFinal planets:");
        for (i, body) in system.bodies.iter().enumerate() {
            let mass_earth = Mass::from_solar_masses(body.mass).to_earth_masses();
            let r = body.position.coords.magnitude();
            let v = body.velocity.magnitude();

            println!(
                "  Planet {}: mass={:.3} M⊕, r={:.3} AU, v={:.3} AU/yr",
                i + 1,
                mass_earth,
                r,
                v
            );
        }
    }

    println!("\n{}", "=".repeat(60));
    println!("Demo complete!");
}
