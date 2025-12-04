//! Generate planets CSV with individual planet data
//!
//! Usage: cargo run -p planetary --example generate_planets
//!
//! Output: planets.csv with one row per planet

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use stellar::solar_analog;

use planetary::from_star;

fn main() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let n_systems = 1000;

    // CSV header
    println!(
        "system_id,planet_idx,mass_earth,radius_earth,sma_au,ecc,inc_deg,eq_temp_k,class,type"
    );

    for system_id in 0..n_systems {
        let star = solar_analog();
        let system = from_star(&mut rng, &star);

        for (planet_idx, planet) in system.planets.iter().enumerate() {
            println!(
                "{},{},{:.4},{:.4},{:.4},{:.4},{:.2},{:.0},{},{}",
                system_id,
                planet_idx,
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

    eprintln!("Generated planets from {} systems", n_systems);
}
