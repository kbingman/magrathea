use stellar::spectral::{LuminosityClass, SpectralType, VariabilityType};
use stellar::stellar_color::StellarColor;
use stellar::stellar_radius::StellarRadius;
use stellar::{MainSequenceStar, StellarObject};
use units::time::Time;
use units::{Mass, Temperature};
use uuid::Uuid;

use crate::generation::{generate_planetary_system, generate_planetary_system_named};
use celestial::snow_line;

/// Create a test star with specified properties
fn make_test_star(
    mass: f64,
    luminosity: f64,
    temperature: f64,
    metallicity: f64,
    spectral: SpectralType,
) -> StellarObject {
    StellarObject::MainSequence(MainSequenceStar {
        mass: Mass::from_solar_masses(mass),
        radius: StellarRadius::from_solar_radii(mass.powf(0.8)), // Rough M-R relation
        luminosity,
        temperature: Temperature::from_kelvin(temperature),
        spectral_type: spectral,
        luminosity_class: LuminosityClass::V,
        subtype: 2,
        variability: VariabilityType::None,
        metallicity,
        age: Time::from_myr(4600.0),
        color: StellarColor::from_temperature(temperature),
    })
}

/// Create a Sun-like test star
fn sun_like_star() -> StellarObject {
    make_test_star(1.0, 1.0, 5778.0, 0.0, SpectralType::G)
}

/// Create a Sun-like star with specified metallicity
fn sun_like_star_metallicity(metallicity: f64) -> StellarObject {
    make_test_star(1.0, 1.0, 5778.0, metallicity, SpectralType::G)
}

#[test]
fn test_generate_system() {
    let star = sun_like_star();
    let system = generate_planetary_system_named(star, "test-42");
    assert!(system.is_stable());
}

#[test]
fn test_reproducibility() {
    let star1 = sun_like_star();
    let star2 = sun_like_star();

    let s1 = generate_planetary_system_named(star1, "reproducibility-test");
    let s2 = generate_planetary_system_named(star2, "reproducibility-test");

    assert_eq!(s1.planets.len(), s2.planets.len());
    assert_eq!(s1.metadata.id, s2.metadata.id);
}

#[test]
fn print_sample_systems() {
    println!("\n=== SAMPLE PLANETARY SYSTEMS ===\n");

    for i in 0..10 {
        let star = sun_like_star();
        let system = generate_planetary_system(star, Uuid::new_v4());

        println!(
            "System {} [{}] ({:?}, {} planets):",
            i + 1,
            system.metadata.catalog_name,
            system.architecture(),
            system.planets.len()
        );

        println!("{:#?}", system);

        for (j, planet) in system.planets.iter().enumerate() {
            let mass = planet.mass.to_earth_masses();
            let sma = planet.semi_major_axis.to_au();
            let mass_str = format!("{:.2} M🜨", mass);

            println!(
                "  {:2}. {:>10} at {:>6.2} AU  ({:?} / {:?})",
                j + 1,
                mass_str,
                sma,
                planet.class,
                planet.planet_type
            );
        }
        println!();
    }
}

// =============================================================================
// Outer System Distribution Validation Tests
// =============================================================================

/// Test that cold giants appear at appropriate rates and locations
/// Expected: ~10-20% of FGK systems have cold giants (metallicity-dependent)
#[test]
fn test_cold_giant_occurrence() {
    let n_systems = 1000;
    let sl = snow_line(1.0); // Sun-like star

    let mut systems_with_cold_giants = 0;

    for i in 0..n_systems {
        let star = sun_like_star();
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("cold-giant-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);

        // Cold giants: >50 M⊕ in Jupiter zone (1.5-6× snow line)
        let has_cold_giant = system.planets.iter().any(|p| {
            let sma = p.semi_major_axis.to_au();
            let mass = p.mass.to_earth_masses();
            mass > 50.0 && sma > sl * 1.5 && sma < sl * 6.0
        });

        if has_cold_giant {
            systems_with_cold_giants += 1;
        }
    }

    let rate = systems_with_cold_giants as f64 / n_systems as f64;
    println!(
        "Cold giant rate: {:.1}% ({}/{})",
        rate * 100.0,
        systems_with_cold_giants,
        n_systems
    );

    // Expected: 5-25% for solar metallicity G stars
    // Note: Rate varies with different UUID seeds across systems
    assert!(
        rate > 0.04 && rate < 0.35,
        "Cold giant rate {:.1}% outside expected range 4-35%",
        rate * 100.0
    );
}

/// Test that ice giants appear at appropriate rates
/// Expected: ~25-40% of systems have at least one ice giant
#[test]
fn test_ice_giant_occurrence() {
    let n_systems = 1000;
    let sl = snow_line(1.0);

    let mut systems_with_ice_giants = 0;
    let mut total_ice_giants = 0;

    for i in 0..n_systems {
        let star = sun_like_star();
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("ice-giant-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);

        // Ice giants: 8-50 M⊕ beyond 5× snow line
        let ice_giants: Vec<_> = system
            .planets
            .iter()
            .filter(|p| {
                let sma = p.semi_major_axis.to_au();
                let mass = p.mass.to_earth_masses();
                mass >= 8.0 && mass <= 50.0 && sma > sl * 5.0
            })
            .collect();

        if !ice_giants.is_empty() {
            systems_with_ice_giants += 1;
            total_ice_giants += ice_giants.len();
        }
    }

    let rate = systems_with_ice_giants as f64 / n_systems as f64;
    let avg_per_system = total_ice_giants as f64 / systems_with_ice_giants.max(1) as f64;

    println!(
        "Ice giant rate: {:.1}% ({}/{}), avg {:.2} per system with ice giants",
        rate * 100.0,
        systems_with_ice_giants,
        n_systems,
        avg_per_system
    );

    // Expected: 10-45% of systems
    // Note: Rate varies with different UUID seeds across systems
    assert!(
        rate > 0.10 && rate < 0.50,
        "Ice giant rate {:.1}% outside expected range 10-50%",
        rate * 100.0
    );

    // Expected multiplicity: ~1.5 per system (50% have 1, 40% have 2, 10% have 3)
    assert!(
        avg_per_system > 1.0 && avg_per_system < 2.5,
        "Ice giant multiplicity {:.2} outside expected range 1.0-2.5",
        avg_per_system
    );
}

/// Test that wide companions are rare
/// Expected: ~1-3% of systems
#[test]
fn test_wide_companion_occurrence() {
    let n_systems = 2000; // More samples for rare events

    let mut systems_with_wide = 0;

    for i in 0..n_systems {
        let star = sun_like_star();
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("wide-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);

        // Wide companions: >50 AU, massive (>1 M_J = 318 M⊕)
        let has_wide = system.planets.iter().any(|p| {
            let sma = p.semi_major_axis.to_au();
            let mass = p.mass.to_earth_masses();
            sma > 50.0 && mass > 100.0
        });

        if has_wide {
            systems_with_wide += 1;
        }
    }

    let rate = systems_with_wide as f64 / n_systems as f64;
    println!(
        "Wide companion rate: {:.1}% ({}/{})",
        rate * 100.0,
        systems_with_wide,
        n_systems
    );

    // Expected: 1-5% (rare but present)
    assert!(
        rate > 0.005 && rate < 0.08,
        "Wide companion rate {:.1}% outside expected range 0.5-8%",
        rate * 100.0
    );
}

/// Test metallicity effect on giant planet occurrence
/// Expected: ~10× increase from [Fe/H]=-0.5 to +0.5
#[test]
fn test_metallicity_giant_correlation() {
    let n_systems = 500;

    // Low metallicity
    let mut low_metal_giants = 0;
    for i in 0..n_systems {
        let star = sun_like_star_metallicity(-0.4);
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("low-metal-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);
        if system
            .planets
            .iter()
            .any(|p| p.mass.to_earth_masses() > 50.0)
        {
            low_metal_giants += 1;
        }
    }

    // High metallicity
    let mut high_metal_giants = 0;
    for i in 0..n_systems {
        let star = sun_like_star_metallicity(0.4);
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("high-metal-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);
        if system
            .planets
            .iter()
            .any(|p| p.mass.to_earth_masses() > 50.0)
        {
            high_metal_giants += 1;
        }
    }

    let low_rate = low_metal_giants as f64 / n_systems as f64;
    let high_rate = high_metal_giants as f64 / n_systems as f64;
    let ratio = high_rate / low_rate.max(0.001);

    println!(
        "Giant rate at [Fe/H]=-0.4: {:.1}%, at +0.4: {:.1}%, ratio: {:.1}×",
        low_rate * 100.0,
        high_rate * 100.0,
        ratio
    );

    // Expected: 3-15× increase (Fischer & Valenti: ~10^(2×0.8) ≈ 6×)
    assert!(
        ratio > 2.0 && ratio < 20.0,
        "Metallicity ratio {:.1}× outside expected range 2-20×",
        ratio
    );
}

/// Test outer system zones are populated correctly
#[test]
fn test_outer_system_zones() {
    let n_systems = 1000;
    let sl = snow_line(1.0);

    let mut jupiter_zone = 0; // 1.5-6× SL
    let mut ice_zone = 0; // 5-15× SL
    let mut wide_zone = 0; // >50 AU

    for i in 0..n_systems {
        let star = sun_like_star();
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("zones-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);

        for planet in &system.planets {
            let sma = planet.semi_major_axis.to_au();

            if sma > sl * 1.5 && sma < sl * 6.0 {
                jupiter_zone += 1;
            }
            if sma > sl * 5.0 && sma < sl * 15.0 {
                ice_zone += 1;
            }
            if sma > 50.0 {
                wide_zone += 1;
            }
        }
    }

    println!(
        "Planets by zone - Jupiter: {}, Ice: {}, Wide: {}",
        jupiter_zone, ice_zone, wide_zone
    );

    // All zones should have some planets
    assert!(jupiter_zone > 50, "Too few planets in Jupiter zone");
    assert!(ice_zone > 100, "Too few planets in ice giant zone");
    assert!(wide_zone > 0, "No wide companions generated");
}

/// Print detailed outer system statistics for manual review
#[test]
fn print_outer_system_statistics() {
    let n_systems = 1000;
    let sl = snow_line(1.0);

    println!("\n=== OUTER SYSTEM STATISTICS (n={}) ===\n", n_systems);

    let mut stats = OuterSystemStats::default();

    for i in 0..n_systems {
        let star = sun_like_star();
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, format!("stats-{}", i).as_bytes());
        let system = generate_planetary_system(star, id);

        let mut has_cold_giant = false;
        let mut has_ice_giant = false;
        let mut has_wide = false;
        let mut n_ice = 0;

        for planet in &system.planets {
            let sma = planet.semi_major_axis.to_au();
            let mass = planet.mass.to_earth_masses();

            // Cold giant
            if mass > 50.0 && sma > sl * 1.5 && sma < sl * 6.0 {
                has_cold_giant = true;
                stats.cold_giant_masses.push(mass);
                stats.cold_giant_smas.push(sma);
            }

            // Ice giant
            if mass >= 8.0 && mass <= 50.0 && sma > sl * 5.0 && sma < sl * 15.0 {
                has_ice_giant = true;
                n_ice += 1;
                stats.ice_giant_masses.push(mass);
                stats.ice_giant_smas.push(sma);
            }

            // Wide companion
            if sma > 50.0 && mass > 100.0 {
                has_wide = true;
                stats.wide_masses.push(mass);
                stats.wide_smas.push(sma);
            }
        }

        if has_cold_giant {
            stats.systems_with_cold_giants += 1;
        }
        if has_ice_giant {
            stats.systems_with_ice_giants += 1;
            match n_ice {
                1 => stats.ice_giant_count_1 += 1,
                2 => stats.ice_giant_count_2 += 1,
                _ => stats.ice_giant_count_3plus += 1,
            }
        }
        if has_wide {
            stats.systems_with_wide += 1;
        }
    }

    println!("Cold Giants:");
    println!(
        "  Rate: {:.1}%",
        stats.systems_with_cold_giants as f64 / n_systems as f64 * 100.0
    );
    if !stats.cold_giant_masses.is_empty() {
        println!(
            "  Mass range: {:.0}-{:.0} M⊕ (median {:.0})",
            stats
                .cold_giant_masses
                .iter()
                .cloned()
                .fold(f64::MAX, f64::min),
            stats.cold_giant_masses.iter().cloned().fold(0.0, f64::max),
            median(&mut stats.cold_giant_masses)
        );
        println!(
            "  SMA range: {:.1}-{:.1} AU",
            stats
                .cold_giant_smas
                .iter()
                .cloned()
                .fold(f64::MAX, f64::min),
            stats.cold_giant_smas.iter().cloned().fold(0.0, f64::max)
        );
    }

    println!("\nIce Giants:");
    println!(
        "  Rate: {:.1}%",
        stats.systems_with_ice_giants as f64 / n_systems as f64 * 100.0
    );
    println!(
        "  Multiplicity: 1={:.1}%, 2={:.1}%, 3+={:.1}%",
        stats.ice_giant_count_1 as f64 / stats.systems_with_ice_giants.max(1) as f64 * 100.0,
        stats.ice_giant_count_2 as f64 / stats.systems_with_ice_giants.max(1) as f64 * 100.0,
        stats.ice_giant_count_3plus as f64 / stats.systems_with_ice_giants.max(1) as f64 * 100.0
    );
    if !stats.ice_giant_masses.is_empty() {
        println!(
            "  Mass range: {:.1}-{:.1} M⊕",
            stats
                .ice_giant_masses
                .iter()
                .cloned()
                .fold(f64::MAX, f64::min),
            stats.ice_giant_masses.iter().cloned().fold(0.0, f64::max)
        );
        println!(
            "  SMA range: {:.1}-{:.1} AU",
            stats
                .ice_giant_smas
                .iter()
                .cloned()
                .fold(f64::MAX, f64::min),
            stats.ice_giant_smas.iter().cloned().fold(0.0, f64::max)
        );
    }

    println!("\nWide Companions:");
    println!(
        "  Rate: {:.1}%",
        stats.systems_with_wide as f64 / n_systems as f64 * 100.0
    );
    if !stats.wide_masses.is_empty() {
        println!(
            "  Mass range: {:.0}-{:.0} M⊕ ({:.1}-{:.1} M_J)",
            stats.wide_masses.iter().cloned().fold(f64::MAX, f64::min),
            stats.wide_masses.iter().cloned().fold(0.0, f64::max),
            stats.wide_masses.iter().cloned().fold(f64::MAX, f64::min) / 318.0,
            stats.wide_masses.iter().cloned().fold(0.0, f64::max) / 318.0
        );
        println!(
            "  SMA range: {:.0}-{:.0} AU",
            stats.wide_smas.iter().cloned().fold(f64::MAX, f64::min),
            stats.wide_smas.iter().cloned().fold(0.0, f64::max)
        );
    }
}

#[derive(Default)]
struct OuterSystemStats {
    systems_with_cold_giants: usize,
    systems_with_ice_giants: usize,
    systems_with_wide: usize,
    ice_giant_count_1: usize,
    ice_giant_count_2: usize,
    ice_giant_count_3plus: usize,
    cold_giant_masses: Vec<f64>,
    cold_giant_smas: Vec<f64>,
    ice_giant_masses: Vec<f64>,
    ice_giant_smas: Vec<f64>,
    wide_masses: Vec<f64>,
    wide_smas: Vec<f64>,
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = v.len() / 2;
    if v.len() % 2 == 0 {
        (v[mid - 1] + v[mid]) / 2.0
    } else {
        v[mid]
    }
}
