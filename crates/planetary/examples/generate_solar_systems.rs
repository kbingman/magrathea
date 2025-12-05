//! Generate planetary systems around solar analog stars and export to CSV
//!
//! Usage: cargo run -p planetary --example generate_solar_systems
//!
//! Output: solar_systems.csv in current directory

use stellar::solar_analog;

use planetary::from_star;

fn main() {
    let n_systems = 1000;

    // CSV header
    println!(
        "system_id,catalog_name,spectral_type,stellar_mass,stellar_luminosity,stellar_temp,metallicity,\
         architecture,n_planets,n_rocky,n_transitional,n_volatile,n_giant,\
         has_hz_planet,innermost_au,outermost_au"
    );

    for i in 0..n_systems {
        let star = solar_analog();
        let system = from_star(&star);

        let n_rocky = system
            .planets
            .iter()
            .filter(|p| p.class == planetary::planet_class::PlanetClass::Rocky)
            .count();
        let n_transitional = system
            .planets
            .iter()
            .filter(|p| p.class == planetary::planet_class::PlanetClass::Transitional)
            .count();
        let n_volatile = system
            .planets
            .iter()
            .filter(|p| p.class == planetary::planet_class::PlanetClass::Volatile)
            .count();
        let n_giant = system
            .planets
            .iter()
            .filter(|p| p.class == planetary::planet_class::PlanetClass::Giant)
            .count();

        let has_hz_planet = !system.habitable_zone_planets().is_empty();

        let innermost = system
            .planets
            .first()
            .map(|p| p.semi_major_axis.to_au())
            .unwrap_or(0.0);
        let outermost = system
            .planets
            .last()
            .map(|p| p.semi_major_axis.to_au())
            .unwrap_or(0.0);

        println!(
            "{},{},{},{:.4},{:.6},{:.0},{:.3},{},{},{},{},{},{},{},{:.4},{:.4}",
            i,
            system.metadata.catalog_name(),
            system.spectral_type(),
            system.effective_mass(),
            system.total_luminosity(),
            system.primary_star().temperature(),
            system.metallicity(),
            system.architecture(),
            system.planets.len(),
            n_rocky,
            n_transitional,
            n_volatile,
            n_giant,
            has_hz_planet,
            innermost,
            outermost,
        );
    }

    eprintln!("Generated {} solar analog systems", n_systems);
}
