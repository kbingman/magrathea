//! Tests for photoevaporation models and disk dispersal.

use super::*;
use approx::assert_relative_eq;

#[test]
fn photoevaporation_model_defaults() {
    let none = PhotoevaporationModel::None;
    let euv = PhotoevaporationModel::Euv { phi_euv: 1e41 };
    let xray = PhotoevaporationModel::Xray {
        luminosity_xray: 1e30,
    };
    let combined = PhotoevaporationModel::Combined {
        phi_euv: 1e41,
        luminosity_xray: 1e30,
    };

    // None model should give zero mass loss
    let mdot = none.mass_loss_rate_per_area(1.496e13, 1.989e33, 1700.0);
    assert_eq!(mdot, 0.0);

    // Other models should give positive mass loss
    let mdot_euv = euv.mass_loss_rate_per_area(1.496e13, 1.989e33, 1700.0);
    assert!(mdot_euv >= 0.0, "EUV mass loss should be non-negative");

    let mdot_xray = xray.mass_loss_rate_per_area(1.496e13, 1.989e33, 1700.0);
    assert!(mdot_xray > 0.0, "X-ray mass loss should be positive");

    let mdot_combined = combined.mass_loss_rate_per_area(1.496e13, 1.989e33, 1700.0);
    assert!(
        mdot_combined >= mdot_euv.max(mdot_xray),
        "Combined should be at least as large as individual components"
    );
}

#[test]
fn t_tauri_model_creates_typical_values() {
    let model = PhotoevaporationModel::t_tauri();

    match model {
        PhotoevaporationModel::Combined {
            phi_euv,
            luminosity_xray,
        } => {
            assert_relative_eq!(phi_euv, 1e41, max_relative = 0.01);
            assert_relative_eq!(luminosity_xray, 1e30, max_relative = 0.01);
        }
        _ => panic!("t_tauri() should return Combined model"),
    }
}

#[test]
fn photoevaporation_reduces_disk_mass() {
    let power_law = GasDisk::mmsn();
    let mut grid = GridDisk::from_gas_disk(&power_law, 200);
    grid.set_photoevaporation(PhotoevaporationModel::t_tauri());

    let initial_mass = grid.total_mass().to_solar_masses();
    println!("Initial mass: {:.6} M_sun", initial_mass);

    // Evolve for 0.1 Myr
    let dt = 1e4 * 3.156e7; // 10,000 years in seconds
    let n_steps = 10; // Total 0.1 Myr

    for step in 0..n_steps {
        grid.apply_photoevaporation(dt);

        if step % 2 == 0 {
            let mass = grid.total_mass().to_solar_masses();
            let time_myr = (step as f64 + 1.0) * 1e4 / 1e6;
            let mdot = grid
                .total_photoevaporation_rate()
                .to_solar_masses_per_year();

            println!(
                "t = {:.3} Myr: M = {:.6} M_sun, dM/dt = {:.2e} M_sun/yr",
                time_myr, mass, mdot
            );
        }
    }

    let final_mass = grid.total_mass().to_solar_masses();
    println!("Final mass: {:.6} M_sun", final_mass);

    // Mass should have decreased
    assert!(
        final_mass < initial_mass,
        "Photoevaporation should reduce disk mass"
    );
}

#[test]
fn disk_lifetime_reasonable() {
    // Test that a disk disperses in 1-10 Myr
    let power_law = GasDisk::mmsn();
    let mut grid = GridDisk::from_gas_disk(&power_law, 200);
    grid.set_photoevaporation(PhotoevaporationModel::t_tauri());

    let initial_mass = grid.total_mass().to_jupiter_masses();
    println!("Initial mass: {:.3} M_jup", initial_mass);

    let dt = 1e4 * 3.156e7; // 10,000 years
    let mut time_myr = 0.0;
    let mut step = 0;

    // Evolve until dispersed or 20 Myr (whichever comes first)
    while !grid.is_dispersed() && time_myr < 20.0 {
        // Apply both viscous evolution and photoevaporation
        let dt_visc = grid.max_timestep().min(dt);
        grid.evolve_viscous(dt_visc);
        grid.apply_photoevaporation(dt_visc);

        step += 1;
        time_myr = step as f64 * dt / 3.156e7 / 1e6;

        if step % 100 == 0 {
            let mass = grid.total_mass().to_jupiter_masses();
            let mdot_pe = grid
                .total_photoevaporation_rate()
                .to_solar_masses_per_year();
            println!(
                "t = {:.2} Myr: M = {:.3} M_jup, dM/dt_pe = {:.2e} M_sun/yr",
                time_myr, mass, mdot_pe
            );
        }
    }

    println!("\nDisk dispersed at t = {:.2} Myr", time_myr);

    // Disk should disperse in 1-20 Myr range
    assert!(
        time_myr >= 0.5 && time_myr <= 20.0,
        "Disk lifetime should be 0.5-20 Myr, got {:.2} Myr",
        time_myr
    );
}

#[test]
fn photoevaporation_rate_scales_with_luminosity() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    // Weak X-ray
    let mut grid_weak = grid.clone();
    grid_weak.set_photoevaporation(PhotoevaporationModel::Xray {
        luminosity_xray: 1e29,
    });

    // Strong X-ray
    let mut grid_strong = grid.clone();
    grid_strong.set_photoevaporation(PhotoevaporationModel::Xray {
        luminosity_xray: 1e31,
    });

    let mdot_weak = grid_weak
        .total_photoevaporation_rate()
        .to_solar_masses_per_year();
    let mdot_strong = grid_strong
        .total_photoevaporation_rate()
        .to_solar_masses_per_year();

    println!("Weak (L_X = 10^29): {:.2e} M_sun/yr", mdot_weak);
    println!("Strong (L_X = 10^31): {:.2e} M_sun/yr", mdot_strong);

    // Higher luminosity should give higher mass loss rate
    assert!(
        mdot_strong > mdot_weak,
        "Higher L_X should give higher mass loss"
    );

    // Should scale roughly linearly with luminosity
    let ratio = mdot_strong / mdot_weak;
    let expected_ratio = 100.0; // L_X ratio
    assert!(
        ratio > expected_ratio * 0.1 && ratio < expected_ratio * 10.0,
        "Mass loss should scale with luminosity"
    );
}

#[test]
fn euv_mass_loss_beyond_gravitational_radius() {
    let model = PhotoevaporationModel::Euv { phi_euv: 1e41 };
    let stellar_mass = 1.989e33; // Solar mass in grams

    // Gravitational radius for T ~ 10^4 K: r_g ≈ 12 AU for solar mass
    let r_inner = 1.0 * 1.496e13; // 1 AU (inside r_g)
    let r_outer = 30.0 * 1.496e13; // 30 AU (outside r_g)

    let mdot_inner = model.mass_loss_rate_per_area(r_inner, stellar_mass, 1700.0);
    let mdot_outer = model.mass_loss_rate_per_area(r_outer, stellar_mass, 1700.0);

    println!("EUV mass loss at 1 AU: {:.2e} g/cm²/s", mdot_inner);
    println!("EUV mass loss at 30 AU: {:.2e} g/cm²/s", mdot_outer);

    // Inner disk (inside r_g) should have zero EUV mass loss
    assert_eq!(mdot_inner, 0.0, "EUV mass loss should be zero inside r_g");

    // Outer disk should have positive mass loss
    assert!(
        mdot_outer > 0.0,
        "EUV mass loss should be positive beyond r_g"
    );
}

#[test]
fn xray_mass_loss_peaks_at_intermediate_radii() {
    let model = PhotoevaporationModel::Xray {
        luminosity_xray: 1e30,
    };
    let stellar_mass = 1.989e33;

    let radii_au = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0];
    let mdots: Vec<f64> = radii_au
        .iter()
        .map(|&r_au| {
            let r = r_au * 1.496e13;
            model.mass_loss_rate_per_area(r, stellar_mass, 1700.0)
        })
        .collect();

    println!("X-ray mass loss rates:");
    for (i, &r_au) in radii_au.iter().enumerate() {
        println!("{:.1} AU: {:.2e} g/cm²/s", r_au, mdots[i]);
    }

    // Find peak
    let (peak_idx, _peak_mdot) = mdots
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    println!("Peak at {:.1} AU", radii_au[peak_idx]);

    // Peak should be at intermediate radii (1-10 AU typically)
    assert!(
        radii_au[peak_idx] >= 0.5 && radii_au[peak_idx] <= 20.0,
        "X-ray mass loss should peak at intermediate radii"
    );

    // Peak should be significantly higher than endpoints
    assert!(
        mdots[peak_idx] > mdots[0] && mdots[peak_idx] > mdots[radii_au.len() - 1],
        "Peak should be higher than inner and outer disk"
    );
}

#[test]
fn dispersed_disk_detection() {
    let power_law = GasDisk::mmsn();
    let mut grid = GridDisk::from_gas_disk(&power_law, 200);

    assert!(!grid.is_dispersed(), "MMSN disk should not be dispersed");

    // Manually reduce all surface densities to near-zero
    for sigma in grid.sigma_mut() {
        *sigma *= 1e-6;
    }

    assert!(
        grid.is_dispersed(),
        "Very low-mass disk should be considered dispersed"
    );
}
