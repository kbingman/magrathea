//! Planetary system generation pipeline

use rand::Rng;
use rand_chacha::ChaChaRng;
use units::{Length, Mass};

use crate::composition::Composition;
use crate::planet::Planet;
use crate::planet_class::PlanetClass;
use crate::sampling::{
    period_to_semi_major_axis, sample_eccentricity, sample_inclination, sample_orbital_period,
    sample_planet_mass,
};
use crate::system::{HabitableZone, PlanetarySystem, SystemArchitecture, snow_line};

/// Generate a complete planetary system
pub fn generate_planetary_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_temperature: f64,
    stellar_metallicity: f64,
    spectral_type: &str,
) -> PlanetarySystem {
    let architecture = SystemArchitecture::sample(rng, spectral_type, stellar_metallicity);

    let planets = match architecture {
        SystemArchitecture::CompactMulti => {
            generate_compact_system(rng, stellar_mass, stellar_luminosity, stellar_metallicity)
        }
        SystemArchitecture::Mixed => {
            generate_mixed_system(rng, stellar_mass, stellar_luminosity, stellar_metallicity)
        }
        SystemArchitecture::GiantDominated => {
            generate_giant_system(rng, stellar_mass, stellar_luminosity, stellar_metallicity)
        }
        SystemArchitecture::Sparse => {
            generate_sparse_system(rng, stellar_mass, stellar_luminosity, stellar_metallicity)
        }
    };

    PlanetarySystem::new(
        stellar_mass,
        stellar_luminosity,
        stellar_temperature,
        stellar_metallicity,
        spectral_type.to_string(),
        planets,
        architecture,
    )
}

fn generate_compact_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_metallicity: f64,
) -> Vec<Planet> {
    let n_planets: usize = rng.random_range(4..=7);
    let hz = HabitableZone::from_luminosity(stellar_luminosity);

    let inner_au = 0.01 * stellar_luminosity.sqrt();
    let outer_au = hz.outer_edge * 1.5;

    generate_spaced_planets(
        rng,
        n_planets,
        inner_au,
        outer_au,
        stellar_mass,
        stellar_luminosity,
        stellar_metallicity,
        0.01..5.0,
    )
}

fn generate_mixed_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_metallicity: f64,
) -> Vec<Planet> {
    let sl = snow_line(stellar_luminosity);

    // Inner terrestrial zone
    let n_inner: usize = rng.random_range(1..=4);
    let inner_planets = generate_spaced_planets(
        rng,
        n_inner,
        0.3 * stellar_luminosity.sqrt(),
        sl * 0.8,
        stellar_mass,
        stellar_luminosity,
        stellar_metallicity,
        0.05..5.0,
    );

    // Outer giant zone
    let giant_prob = 0.3 * 10.0_f64.powf(stellar_metallicity);
    let outer_planets = if rng.random::<f64>() < giant_prob {
        let n_outer: usize = rng.random_range(1..=2);
        generate_spaced_planets(
            rng,
            n_outer,
            sl * 1.5,
            sl * 10.0,
            stellar_mass,
            stellar_luminosity,
            stellar_metallicity,
            50.0..500.0,
        )
    } else {
        vec![]
    };

    [inner_planets, outer_planets].concat()
}

fn generate_giant_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_metallicity: f64,
) -> Vec<Planet> {
    let sl = snow_line(stellar_luminosity);
    let hot_jupiter = rng.random::<f64>() < 0.2;

    let mut planets = if hot_jupiter {
        let mass = sample_giant_mass(rng);
        let period = rng.random_range(2.0..10.0);
        let sma = period_to_semi_major_axis(period, stellar_mass);
        vec![create_planet(
            rng,
            mass,
            sma,
            stellar_mass,
            stellar_luminosity,
        )]
    } else {
        let n_giants: usize = rng.random_range(1..=2);
        generate_spaced_planets(
            rng,
            n_giants,
            sl * 1.5,
            sl * 15.0,
            stellar_mass,
            stellar_luminosity,
            stellar_metallicity,
            100.0..1000.0,
        )
    };

    // Maybe surviving inner planets
    if !hot_jupiter && rng.random::<f64>() < 0.3 {
        let sl = snow_line(stellar_luminosity);
        let n_inner = rng.random_range(1..=2);
        let inner = generate_spaced_planets(
            rng,
            n_inner,
            0.3 * stellar_luminosity.sqrt(),
            sl * 0.5,
            stellar_mass,
            stellar_luminosity,
            stellar_metallicity,
            0.1..3.0,
        );
        planets = [inner, planets].concat();
    }

    planets
}

fn generate_sparse_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_metallicity: f64,
) -> Vec<Planet> {
    if rng.random::<f64>() < 0.5 {
        return vec![];
    }

    let mass = sample_planet_mass(rng, stellar_metallicity);
    let class = PlanetClass::from_earth_masses(mass);
    let period = sample_orbital_period(rng, &class);
    let sma = period_to_semi_major_axis(period, stellar_mass);

    vec![create_planet(
        rng,
        mass,
        sma,
        stellar_mass,
        stellar_luminosity,
    )]
}

fn generate_spaced_planets(
    rng: &mut ChaChaRng,
    n_planets: usize,
    inner_au: f64,
    outer_au: f64,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_metallicity: f64,
    mass_range: std::ops::Range<f64>,
) -> Vec<Planet> {
    if n_planets == 0 {
        return vec![];
    }

    let mut last_attempt = Vec::new();

    for _ in 0..50 {
        let mut planets = Vec::with_capacity(n_planets);

        let log_inner = inner_au.ln();
        let log_outer = outer_au.ln();
        let log_spacing = (log_outer - log_inner) / n_planets as f64;

        for i in 0..n_planets {
            let log_base = log_inner + log_spacing * (i as f64 + 0.5);
            let scatter = rng.random_range(-0.3..0.3) * log_spacing;
            let sma = (log_base + scatter).exp();

            let mass =
                sample_mass_in_range(rng, mass_range.start, mass_range.end, stellar_metallicity);
            planets.push(create_planet(
                rng,
                mass,
                sma,
                stellar_mass,
                stellar_luminosity,
            ));
        }

        planets.sort_by(|a, b| {
            a.semi_major_axis
                .to_au()
                .partial_cmp(&b.semi_major_axis.to_au())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if is_stable(&planets, stellar_mass) {
            return planets;
        }

        last_attempt = planets;
    }

    last_attempt
}

/// Sample mass with metallicity-dependent weighting
///
/// Higher metallicity shifts the distribution toward higher masses,
/// reflecting the giant planet-metallicity correlation (Fischer & Valenti 2005).
/// Base exponent 0.7 biases toward lower masses; metallicity adjusts this.
fn sample_mass_in_range(rng: &mut ChaChaRng, min: f64, max: f64, metallicity: f64) -> f64 {
    let log_min = min.max(0.01).ln();
    let log_max = max.ln();

    // Metallicity effect: higher [Fe/H] → flatter distribution (more high-mass planets)
    // At [Fe/H] = 0 (solar): exponent = 0.7 (baseline)
    // At [Fe/H] = +0.3: exponent ≈ 0.55 (more giants)
    // At [Fe/H] = -0.3: exponent ≈ 0.85 (fewer giants)
    let exponent = (0.7 - metallicity * 0.5).clamp(0.4, 0.95);

    let u: f64 = rng.random();
    (log_min + u.powf(exponent) * (log_max - log_min)).exp()
}

fn sample_giant_mass(rng: &mut ChaChaRng) -> f64 {
    let log_jupiter = 318.0_f64.ln();
    let z: f64 = rng.random_range(-0.5..0.5);
    (log_jupiter + z).exp().clamp(50.0, 2000.0)
}

fn create_planet(
    rng: &mut ChaChaRng,
    mass_earth: f64,
    sma_au: f64,
    stellar_mass: f64,
    stellar_luminosity: f64,
) -> Planet {
    let class = PlanetClass::from_earth_masses(mass_earth);
    let period = (sma_au.powi(3) / stellar_mass).sqrt() * 365.25;
    let eccentricity = sample_eccentricity(rng, &class, period);
    let inclination = sample_inclination(rng, true);
    let composition = sample_composition(rng, &class, sma_au, stellar_luminosity);

    Planet::from_mass(
        Mass::from_earth_masses(mass_earth),
        Length::from_au(sma_au),
        eccentricity,
        inclination,
        composition,
        stellar_luminosity,
        rng,
    )
}

fn sample_composition(
    rng: &mut ChaChaRng,
    class: &PlanetClass,
    sma_au: f64,
    stellar_luminosity: f64,
) -> Composition {
    let sl = snow_line(stellar_luminosity);
    let beyond_snow_line = sma_au > sl;

    match class {
        PlanetClass::Rocky => {
            if beyond_snow_line {
                Composition::new(
                    0.15 + rng.random::<f64>() * 0.10,
                    0.35 + rng.random::<f64>() * 0.10,
                    0.40 + rng.random::<f64>() * 0.20,
                    0.0,
                )
            } else {
                Composition::new(
                    0.25 + rng.random::<f64>() * 0.15,
                    0.60 + rng.random::<f64>() * 0.15,
                    rng.random::<f64>() * 0.05,
                    0.0,
                )
            }
        }
        PlanetClass::Transitional => {
            let water = if beyond_snow_line {
                0.20 + rng.random::<f64>() * 0.30
            } else {
                rng.random::<f64>() * 0.10
            };
            let envelope = rng.random::<f64>() * 0.15;
            Composition::new(
                0.15 + rng.random::<f64>() * 0.10,
                0.40 + rng.random::<f64>() * 0.10,
                water,
                envelope,
            )
        }
        PlanetClass::Volatile => Composition::sample_ice_giant(rng),
        PlanetClass::Giant => Composition::sample_gas_giant(rng),
    }
}

fn is_stable(planets: &[Planet], stellar_mass: f64) -> bool {
    if planets.len() < 2 {
        return true;
    }

    for window in planets.windows(2) {
        let m1 = window[0].mass.to_solar_masses();
        let m2 = window[1].mass.to_solar_masses();
        let a1 = window[0].semi_major_axis.to_au();
        let a2 = window[1].semi_major_axis.to_au();

        let mutual_hill = ((m1 + m2) / (3.0 * stellar_mass)).powf(1.0 / 3.0) * ((a1 + a2) / 2.0);
        if (a2 - a1) / mutual_hill < 8.0 {
            return false;
        }
    }

    true
}
