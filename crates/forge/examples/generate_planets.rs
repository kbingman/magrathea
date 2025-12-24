//! Generate planets CSV with individual planet data
//!
//! Usage: cargo run -p forge --example generate_planets
//!
//! Output: planets.csv with one row per planet

use stellar::solar_analog;

use forge::from_star;

fn main() {
    let n_systems = 1000;

    // CSV header
    println!(
        "system_id,catalog_name,planet_idx,mass_earth,radius_earth,sma_au,ecc,inc_deg,eq_temp_k,class,type"
    );

    for system_id in 0..n_systems {
        let star = solar_analog();
        let system = from_star(&star);

        for (planet_idx, planet) in system.planets.iter().enumerate() {
            println!(
                "{},{},{},{:.4},{:.4},{:.4},{:.4},{:.2},{:.0},{},{}",
                system_id,
                system.metadata.catalog_name,
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
