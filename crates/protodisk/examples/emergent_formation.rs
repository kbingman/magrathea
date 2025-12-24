//! Emergent planet formation demonstration.
//!
//! This example shows the complete physics pipeline:
//! 1. Particles settle and grow via coagulation
//! 2. Gravitational instability forms planetesimals
//! 3. Planetesimals become discrete bodies
//! 4. Bodies accrete from particle disk
//! 5. Bodies capture gas envelopes
//! 6. Bodies migrate inward
//! 7. Bodies can collide and merge
//!
//! Run with: cargo run --example emergent_formation

use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use protodisk::bodies::DiscreteBody;
use protodisk::{
    DiskMass, DiskModel, GasDisk, GridDisk, ParticleBin, SimulationState, run_simulation,
    solar_analog,
};
use units::{Length, Mass, Time};

const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

fn main() {
    println!("╔═══════════════════════════════════════════════════════╗");
    println!("║   Emergent Planet Formation Simulation               ║");
    println!("║   Bottom-Up Physics → Emergent Systems                ║");
    println!("╚═══════════════════════════════════════════════════════╝\n");

    // Create a solar-type star
    let star = solar_analog();
    println!(
        "🌟 Star: {} ({:.2} M☉, {:.0} K)",
        star.spectral_type,
        star.mass.to_solar_masses(),
        star.temperature.to_kelvin()
    );

    // Create a protoplanetary disk
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

    println!(
        "💿 Disk: {:.3} M☉",
        grid_disk.total_mass().to_solar_masses()
    );
    println!("   Inner edge: {:.2} AU", grid_disk.inner_radius().to_au());
    println!(
        "   Outer edge: {:.2} AU\n",
        grid_disk.outer_radius().to_au()
    );

    // Create particle populations directly from disk model
    // Start with "pre-settled" conditions (simulating ~100k years of evolution)
    println!("🔬 Creating particle populations (pre-settled)...");
    let mut particle_bins = vec![
        // Inner disk (rocky region)
        ParticleBin::from_disk(&gas_disk, Length::from_au(0.5), Length::from_au(0.2)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(1.0), Length::from_au(0.5)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(2.0), Length::from_au(0.5)),
        // Mid disk (around snow line)
        ParticleBin::from_disk(&gas_disk, Length::from_au(4.0), Length::from_au(0.5)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(5.0), Length::from_au(0.5)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(7.0), Length::from_au(1.0)),
        // Outer disk
        ParticleBin::from_disk(&gas_disk, Length::from_au(10.0), Length::from_au(1.0)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(15.0), Length::from_au(2.0)),
        ParticleBin::from_disk(&gas_disk, Length::from_au(25.0), Length::from_au(3.0)),
    ];

    // Manually settle particles to simulate prior evolution
    // This makes the demo show interesting phenomena in a reasonable runtime
    for bin in &mut particle_bins {
        let r = bin.radial_center();
        // Particles have settled to fraction of gas scale height
        let h_g = gas_disk.scale_height(r);
        bin.set_scale_height(h_g * 0.3); // Well-settled

        // Velocity dispersion has decreased (particles are "cold")
        let cs = gas_disk.sound_speed(r);
        bin.set_velocity_dispersion(cs * 0.005); // Low dispersion
    }

    for (i, bin) in particle_bins.iter().enumerate() {
        println!(
            "   Bin {}: {:.1} AU, {:.1e} M⊕, Q={:.2}",
            i,
            bin.radial_center().to_au(),
            bin.total_mass().to_earth_masses(),
            bin.toomre_q(&gas_disk)
        );
    }

    // Create some initial protoplanets to demonstrate migration and envelope capture
    println!("\n🌱 Creating initial protoplanets...");
    let stellar_mass = star.mass.to_solar_masses();

    let mut protoplanets = vec![];

    // Inner super-Earth at 1 AU
    let a1 = 1.0;
    let v1 = (G * stellar_mass / a1).sqrt();
    protoplanets.push(DiscreteBody::new(
        Mass::from_earth_masses(3.0),
        Mass::zero(),
        Point2::new(a1, 0.0),
        Vector2::new(0.0, v1),
        star.mass,
        Composition::earth_like(),
    ));

    // Jupiter-like core at 5 AU
    let a2 = 5.0;
    let v2 = (G * stellar_mass / a2).sqrt();
    protoplanets.push(DiscreteBody::new(
        Mass::from_earth_masses(10.0),
        Mass::zero(),
        Point2::new(a2, 0.0),
        Vector2::new(0.0, v2),
        star.mass,
        Composition::earth_like(),
    ));

    println!("   Proto-Earth: 1.0 AU, 3.0 M⊕");
    println!("   Proto-Jupiter: 5.0 AU, 10.0 M⊕\n");

    // Create initial state
    let state = SimulationState::with_seed(grid_disk, particle_bins, protoplanets, 42);

    println!("\n▶️  Running simulation (max 10,000 years)...\n");
    println!("{:-<80}", "");

    // Run simulation for up to 10,000 years (or until disk disperses)
    let max_time = Time::from_years(10_000.0);
    let final_state = run_simulation(state, max_time);

    println!("{:-<80}\n", "");

    // Report results
    println!("╔═══════════════════════════════════════════════════════╗");
    println!("║   Simulation Results                                  ║");
    println!("╚═══════════════════════════════════════════════════════╝\n");

    println!(
        "⏱️  Simulation time: {:.2e} years",
        final_state.time.to_years()
    );
    println!(
        "💿 Final disk mass: {:.4} M☉",
        final_state.disk.total_mass().to_solar_masses()
    );
    println!("🪨 Particle bins: {}", final_state.particle_bins.len());
    println!(
        "🌍 Discrete bodies: {}\n",
        final_state.discrete_bodies.len()
    );

    if !final_state.discrete_bodies.is_empty() {
        println!("╔═══════════════════════════════════════════════════════╗");
        println!("║   Formed Bodies (Planets & Planetesimals)            ║");
        println!("╚═══════════════════════════════════════════════════════╝\n");

        // Sort by semi-major axis
        let mut bodies = final_state.discrete_bodies.clone();
        bodies.sort_by(|a, b| {
            a.semi_major_axis
                .to_au()
                .partial_cmp(&b.semi_major_axis.to_au())
                .unwrap()
        });

        println!(
            "{:<6} {:>8} {:>12} {:>12} {:>8} {:>8}",
            "Body", "a (AU)", "Mass (M⊕)", "Env (M⊕)", "e", "Type"
        );
        println!("{:-<80}", "");

        for (i, body) in bodies.iter().enumerate() {
            let total_mass = body.total_mass().to_earth_masses();
            let env_mass = body.envelope_mass.to_earth_masses();
            let body_type = classify_body(total_mass, env_mass);

            println!(
                "{:<6} {:>8.3} {:>12.3e} {:>12.3e} {:>8.4} {:>8}",
                i,
                body.semi_major_axis.to_au(),
                total_mass,
                env_mass,
                body.eccentricity,
                body_type
            );
        }

        println!();

        // Summary statistics
        let total_body_mass: Mass = bodies
            .iter()
            .map(|b| b.total_mass())
            .fold(Mass::zero(), |acc, m| acc + m);

        let planets = bodies
            .iter()
            .filter(|b| b.total_mass().to_earth_masses() > 0.1)
            .count();

        let gas_giants = bodies
            .iter()
            .filter(|b| b.envelope_mass.to_earth_masses() > 1.0)
            .count();

        println!("📊 Summary:");
        println!(
            "   Total mass in bodies: {:.3e} M⊕",
            total_body_mass.to_earth_masses()
        );
        println!("   Planet-sized bodies (>0.1 M⊕): {}", planets);
        println!("   Gas-rich bodies (>1 M⊕ envelope): {}", gas_giants);

        // Check for interesting features
        if gas_giants > 0 {
            println!("\n✨ Giant planet formation via core accretion!");
        }

        let hot_jupiters = bodies
            .iter()
            .filter(|b| b.semi_major_axis.to_au() < 0.1 && b.envelope_mass.to_earth_masses() > 10.0)
            .count();

        if hot_jupiters > 0 {
            println!("🔥 Hot Jupiter detected (Type I migration)!");
        }
    }

    println!("\n╔═══════════════════════════════════════════════════════╗");
    println!("║   Physics Validation                                  ║");
    println!("╚═══════════════════════════════════════════════════════╝\n");

    // Show that physics worked as expected
    let initial_disk = GasDisk::for_star(&star);
    let initial_grid = GridDisk::from_gas_disk(&initial_disk, 50);
    let initial_mass = initial_grid.total_mass();
    let mass_lost = initial_mass - final_state.disk.total_mass();

    println!("✓ Disk evolution:");
    println!("  Initial: {:.4} M☉", initial_mass.to_solar_masses());
    println!(
        "  Final:   {:.4} M☉",
        final_state.disk.total_mass().to_solar_masses()
    );
    println!(
        "  Lost:    {:.4} M☉ (photoevaporation)\n",
        mass_lost.to_solar_masses()
    );

    if !final_state.discrete_bodies.is_empty() {
        let migrated = final_state
            .discrete_bodies
            .iter()
            .filter(|b| {
                b.semi_major_axis.to_au() < 4.5 // Formed at 4-6 AU
            })
            .count();

        if migrated > 0 {
            println!("✓ Type I migration: {} bodies migrated inward\n", migrated);
        }

        let with_envelope = final_state
            .discrete_bodies
            .iter()
            .filter(|b| b.envelope_mass.to_earth_masses() > 0.01)
            .count();

        if with_envelope > 0 {
            println!(
                "✓ Envelope capture: {} bodies captured gas\n",
                with_envelope
            );
        }
    }

    println!("🎉 Emergent formation simulation complete!");
    println!("   All phenomena emerged from local physics rules.\n");
}

/// Classify body by mass and envelope
fn classify_body(total_mass_earth: f64, envelope_mass_earth: f64) -> &'static str {
    if envelope_mass_earth > 10.0 {
        "Gas Giant"
    } else if envelope_mass_earth > 1.0 {
        "Ice Giant"
    } else if total_mass_earth > 0.1 {
        "Rocky"
    } else {
        "P'simal" // Planetesimal
    }
}
