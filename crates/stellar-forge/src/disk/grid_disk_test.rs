use approx::assert_relative_eq;

use crate::disk::{DiskMass, DiskModel, GasDisk, GridDisk};
use units::Length;

#[test]
fn from_gas_disk_matches_surface_density() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    for r_au in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0] {
        let r = Length::from_au(r_au);
        let sigma_pl = power_law.surface_density(r).to_grams_per_cm2();
        let sigma_grid = grid.surface_density(r).to_grams_per_cm2();

        assert_relative_eq!(sigma_grid, sigma_pl, max_relative = 0.01);
    }
}

#[test]
fn from_gas_disk_matches_temperature() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    for r_au in [0.5, 1.0, 5.0, 10.0] {
        let r = Length::from_au(r_au);
        let t_pl = power_law.temperature(r).to_kelvin();
        let t_grid = grid.temperature(r).to_kelvin();

        assert_relative_eq!(t_grid, t_pl, max_relative = 1e-10);
    }
}

#[test]
fn derived_quantities_match() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let r = Length::from_au(1.0);

    // Orbital frequency
    let omega_pl = power_law.orbital_frequency(r).to_rad_per_sec();
    let omega_grid = grid.orbital_frequency(r).to_rad_per_sec();
    assert_relative_eq!(omega_grid, omega_pl, max_relative = 1e-10);

    // Sound speed
    let cs_pl = power_law.sound_speed(r).to_cm_per_sec();
    let cs_grid = grid.sound_speed(r).to_cm_per_sec();
    assert_relative_eq!(cs_grid, cs_pl, max_relative = 1e-10);

    // Scale height
    let h_pl = power_law.scale_height(r).to_cm();
    let h_grid = grid.scale_height(r).to_cm();
    assert_relative_eq!(h_grid, h_pl, max_relative = 1e-10);

    // Aspect ratio
    let hr_pl = power_law.aspect_ratio(r);
    let hr_grid = grid.aspect_ratio(r);
    assert_relative_eq!(hr_grid, hr_pl, max_relative = 1e-10);
}

#[test]
fn pressure_gradient_numerical_vs_analytical() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let r = Length::from_au(1.0);

    // GasDisk uses analytical: -(p + (3+q)/2) = -2.75 for MMSN
    let d_ln_p_analytical = power_law.pressure_gradient_log(r);

    // GridDisk uses numerical central difference
    let d_ln_p_numerical = grid.pressure_gradient_log(r);

    // Should match within 1%
    assert_relative_eq!(d_ln_p_numerical, d_ln_p_analytical, max_relative = 0.01);
}

#[test]
fn eta_matches() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let r = Length::from_au(1.0);

    let eta_pl = power_law.pressure_gradient_parameter(r);
    let eta_grid = grid.pressure_gradient_parameter(r);

    // Should match within 1% (limited by numerical differentiation)
    assert_relative_eq!(eta_grid, eta_pl, max_relative = 0.01);
}

#[test]
fn total_mass_matches() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let mass_pl = power_law.total_mass().to_solar_masses();
    let mass_grid = grid.total_mass().to_solar_masses();

    // Should match within 1%
    assert_relative_eq!(mass_grid, mass_pl, max_relative = 0.01);
}

#[test]
fn grid_has_correct_bounds() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    assert_relative_eq!(
        grid.inner_radius().to_au(),
        power_law.inner_radius.to_au(),
        max_relative = 1e-10
    );
    assert_relative_eq!(
        grid.outer_radius().to_au(),
        power_law.outer_radius.to_au(),
        max_relative = 1e-10
    );
}

#[test]
fn n_radii_correct() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    assert_eq!(grid.n_radii(), 200);
}

#[test]
fn interpolation_at_grid_points_exact() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 100);

    // At grid points, interpolation should be exact
    for (i, &r) in grid.radii().iter().enumerate() {
        let sigma_interp = grid.surface_density(Length::from_cm(r)).to_grams_per_cm2();
        let sigma_stored = grid.sigma()[i];

        assert_relative_eq!(sigma_interp, sigma_stored, max_relative = 1e-10);
    }
}

// =============================================================================
// Viscous evolution tests
// =============================================================================

#[test]
fn max_timestep_is_positive() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let dt = grid.max_timestep();
    assert!(dt > 0.0, "max_timestep should be positive");
    assert!(dt < 1e20, "max_timestep should be finite");
}

#[test]
fn viscous_evolution_depletes_inner_disk() {
    let power_law = GasDisk::mmsn();
    let mut grid = GridDisk::from_gas_disk(&power_law, 200);

    // Record initial surface density at inner region
    let sigma_inner_initial = grid
        .surface_density(Length::from_au(0.5))
        .to_grams_per_cm2();

    // Evolve for 0.1 viscous times at 1 AU
    let t_visc = grid.viscous_timescale(Length::from_au(1.0));
    let evolution_time = 0.1 * t_visc.to_seconds();

    let dt = grid.max_timestep();
    let n_steps = (evolution_time / dt).ceil() as usize;

    for _ in 0..n_steps {
        grid.evolve_viscous(dt);
    }

    // Inner disk should have depleted
    let sigma_inner_final = grid
        .surface_density(Length::from_au(0.5))
        .to_grams_per_cm2();

    assert!(
        sigma_inner_final < sigma_inner_initial,
        "Inner disk should deplete: initial {:.2e}, final {:.2e}",
        sigma_inner_initial,
        sigma_inner_final
    );
}

#[test]
fn viscous_evolution_spreads_outer_disk() {
    let power_law = GasDisk::mmsn();
    let mut grid = GridDisk::from_gas_disk(&power_law, 200);

    // Record initial surface density at outer region
    let sigma_outer_initial = grid
        .surface_density(Length::from_au(50.0))
        .to_grams_per_cm2();

    // Evolve for 0.1 viscous times at 1 AU
    let t_visc = grid.viscous_timescale(Length::from_au(1.0));
    let evolution_time = 0.1 * t_visc.to_seconds();

    let dt = grid.max_timestep();
    let n_steps = (evolution_time / dt).ceil() as usize;

    for _ in 0..n_steps {
        grid.evolve_viscous(dt);
    }

    // Outer disk should have gained mass (spreading)
    let sigma_outer_final = grid
        .surface_density(Length::from_au(50.0))
        .to_grams_per_cm2();

    assert!(
        sigma_outer_final > sigma_outer_initial,
        "Outer disk should spread: initial {:.2e}, final {:.2e}",
        sigma_outer_initial,
        sigma_outer_final
    );
}

#[test]
fn viscous_timescale_reasonable() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);

    let t_visc = grid.viscous_timescale(Length::from_au(1.0));

    // Should be ~10^5 to 10^6 years for α = 10^(-3)
    let t_yr = t_visc.to_years();
    assert!(
        t_yr > 1e4 && t_yr < 1e7,
        "t_visc at 1 AU: expected ~10^5-10^6 yr, got {:.2e} yr",
        t_yr
    );
}
