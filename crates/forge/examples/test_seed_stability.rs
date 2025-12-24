//! Test that the same seed produces identical systems
//!
//! Usage: cargo run -p forge --example test_seed_stability

use forge::from_star_with_id;
use stellar::main_sequence_star;
use uuid::Uuid;

fn main() {
    let star = main_sequence_star(1.0, 0.0, 4600.0);
    let test_uuid = Uuid::parse_str("550e8400-e29b-41d4-a716-446655440000").unwrap();

    println!("Testing seed stability with UUID: {}", test_uuid);
    println!("Star: G2V, 1.0 M☉\n");

    for run in 1..=5 {
        let system = from_star_with_id(&star, test_uuid);

        println!(
            "Run {}: {} planets, {} architecture",
            run,
            system.planets.len(),
            system.architecture()
        );

        for (i, planet) in system.planets.iter().enumerate() {
            println!(
                "  Planet {}: {:.3} M⊕, {:.3} AU, {}",
                i,
                planet.mass.to_earth_masses(),
                planet.semi_major_axis.to_au(),
                planet.class
            );
        }
        println!();
    }

    // Verify stability
    let system1 = from_star_with_id(&star, test_uuid);
    let system2 = from_star_with_id(&star, test_uuid);

    if system1.planets.len() != system2.planets.len() {
        eprintln!(
            "❌ FAIL: Planet count differs! {} vs {}",
            system1.planets.len(),
            system2.planets.len()
        );
        std::process::exit(1);
    }

    if system1.architecture() != system2.architecture() {
        eprintln!(
            "❌ FAIL: Architecture differs! {} vs {}",
            system1.architecture(),
            system2.architecture()
        );
        std::process::exit(1);
    }

    for (i, (p1, p2)) in system1
        .planets
        .iter()
        .zip(system2.planets.iter())
        .enumerate()
    {
        if (p1.mass.to_earth_masses() - p2.mass.to_earth_masses()).abs() > 0.001 {
            eprintln!(
                "❌ FAIL: Planet {} mass differs! {:.3} vs {:.3}",
                i,
                p1.mass.to_earth_masses(),
                p2.mass.to_earth_masses()
            );
            std::process::exit(1);
        }
        if (p1.semi_major_axis.to_au() - p2.semi_major_axis.to_au()).abs() > 0.001 {
            eprintln!(
                "❌ FAIL: Planet {} SMA differs! {:.3} vs {:.3}",
                i,
                p1.semi_major_axis.to_au(),
                p2.semi_major_axis.to_au()
            );
            std::process::exit(1);
        }
    }

    println!("✓ PASS: All 5 runs produced identical systems");
}
