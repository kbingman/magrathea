//! Generate planets CSV with individual planet data from IMF-sampled stars
//!
//! Usage: cargo run -p forge --example generate_planets_imf
//!
//! Output: planets_imf.csv with one row per planet

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use stellar::sample_main_sequence_star;

use forge::from_star;

fn main() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let n_systems = 1000;

    // CSV header
    println!(
        "system_id,catalog_name,planet_idx,spectral_type,stellar_mass,mass_earth,radius_earth,sma_au,ecc,inc_deg,eq_temp_k,class,type"
    );

    for system_id in 0..n_systems {
        let star = sample_main_sequence_star(&mut rng);
        let system = from_star(&star);

        for (planet_idx, planet) in system.planets.iter().enumerate() {
            println!(
                "{},{},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.2},{:.0},{},{}",
                system_id,
                system.metadata.catalog_name,
                planet_idx,
                system.spectral_type(),
                system.effective_mass(),
                planet.mass.to_earth_masses(),
                planet.radius.to_earth_radii(),
                planet.semi_major_axis.to_au(),
                planet.eccentricity,
                planet.inclination.to_degrees(),
                planet.equilibrium_temp,
                planet.class,
                planet.planet_type,
            );
        }
    }

    eprintln!("Generated planets from {} IMF-sampled systems", n_systems);
}
