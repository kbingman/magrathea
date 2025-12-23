use approx::assert_relative_eq;

use crate::disk::{DiskMass, DiskModel, GasDisk};
use units::Length;

#[test]
fn mmsn_surface_density_at_1au() {
    let disk = GasDisk::mmsn();
    let sigma = disk.surface_density(Length::from_au(1.0));

    // Should be 1700 g/cm² at 1 AU by construction
    assert_relative_eq!(sigma.to_grams_per_cm2(), 1700.0, max_relative = 1e-10);
}

#[test]
fn mmsn_surface_density_at_5au() {
    let disk = GasDisk::mmsn();
    let sigma = disk.surface_density(Length::from_au(5.0));

    // Σ(5 AU) = 1700 × 5^(-1) = 340 g/cm²
    let expected = 1700.0 / 5.0;
    assert_relative_eq!(sigma.to_grams_per_cm2(), expected, max_relative = 1e-10);
}

#[test]
fn mmsn_temperature_at_1au() {
    let disk = GasDisk::mmsn();
    let t = disk.temperature(Length::from_au(1.0));

    // Should be 280 K at 1 AU by construction
    assert_relative_eq!(t.to_kelvin(), 280.0, max_relative = 1e-10);
}

#[test]
fn mmsn_temperature_at_5au() {
    let disk = GasDisk::mmsn();
    let t = disk.temperature(Length::from_au(5.0));

    // T(5 AU) = 280 × 5^(-0.5) ≈ 125 K
    let expected = 280.0 / 5.0_f64.sqrt();
    assert_relative_eq!(t.to_kelvin(), expected, max_relative = 1e-10);
}

#[test]
fn orbital_period_at_1au() {
    let disk = GasDisk::mmsn();
    let period = disk.orbital_period(Length::from_au(1.0));

    // Should be ~1 year
    assert_relative_eq!(period.to_years(), 1.0, max_relative = 0.01);
}

#[test]
fn orbital_period_at_5au() {
    let disk = GasDisk::mmsn();
    let period = disk.orbital_period(Length::from_au(5.0));

    // Kepler's third law: P ∝ a^(3/2), so P(5 AU) ≈ 11.2 yr
    let expected = 5.0_f64.powf(1.5);
    assert_relative_eq!(period.to_years(), expected, max_relative = 0.01);
}

#[test]
fn aspect_ratio_reasonable() {
    let disk = GasDisk::mmsn();

    // h/r should be ~0.03-0.05 at 1 AU, increasing outward for flared disk
    let h_r_1au = disk.aspect_ratio(Length::from_au(1.0));
    let h_r_5au = disk.aspect_ratio(Length::from_au(5.0));

    assert!(
        h_r_1au > 0.02 && h_r_1au < 0.07,
        "h/r at 1 AU: expected 0.03-0.05, got {:.4}",
        h_r_1au
    );

    // For q = 0.5, h/r ∝ r^(0.25), so h/r increases outward
    assert!(
        h_r_5au > h_r_1au,
        "Disk should be flared: h/r at 5 AU ({:.4}) should exceed 1 AU ({:.4})",
        h_r_5au,
        h_r_1au
    );
}

#[test]
fn pressure_gradient_parameter_reasonable() {
    let disk = GasDisk::mmsn();

    // η should be ~0.002-0.005 for typical disk parameters
    let eta_1au = disk.pressure_gradient_parameter(Length::from_au(1.0));
    let eta_5au = disk.pressure_gradient_parameter(Length::from_au(5.0));

    assert!(
        eta_1au > 0.001 && eta_1au < 0.01,
        "η at 1 AU: expected ~0.002-0.005, got {:.5}",
        eta_1au
    );

    // η ∝ (h/r)² ∝ r^(0.5) for q = 0.5, so η increases outward
    assert!(
        eta_5au > eta_1au,
        "η should increase outward: 5 AU ({:.5}) > 1 AU ({:.5})",
        eta_5au,
        eta_1au
    );

    println!("η at 1 AU: {:.5}", eta_1au);
    println!("η at 5 AU: {:.5}", eta_5au);
}

#[test]
fn sub_keplerian_velocity_reasonable() {
    let disk = GasDisk::mmsn();

    // Δv = η × v_K, should be ~50-100 m/s at 1 AU
    let delta_v = disk.sub_keplerian_velocity(Length::from_au(1.0));
    let delta_v_m_s = delta_v.to_cm_per_sec() / 100.0; // convert to m/s

    assert!(
        delta_v_m_s > 30.0 && delta_v_m_s < 150.0,
        "Δv at 1 AU: expected ~50-100 m/s, got {:.1} m/s",
        delta_v_m_s
    );

    println!("Sub-Keplerian velocity at 1 AU: {:.1} m/s", delta_v_m_s);
}

#[test]
fn midplane_density_reasonable() {
    let disk = GasDisk::mmsn();

    // ρ at 1 AU should be ~10^(-10) to 10^(-9) g/cm³
    let rho = disk.midplane_density(Length::from_au(1.0));

    assert!(
        rho.to_grams_per_cm3() > 1e-11 && rho.to_grams_per_cm3() < 1e-8,
        "ρ at 1 AU: expected ~10^(-10) g/cm³, got {:.2e}",
        rho.to_grams_per_cm3()
    );

    println!(
        "Midplane density at 1 AU: {:.2e} g/cm³",
        rho.to_grams_per_cm3()
    );
}

#[test]
fn total_mass_reasonable() {
    let disk = GasDisk::mmsn();
    let mass = disk.total_mass();

    // MMSN total mass should be ~0.01-0.1 solar masses
    let mass_solar = mass.to_solar_masses();

    assert!(
        mass_solar > 0.005 && mass_solar < 0.5,
        "Total disk mass: expected ~0.01-0.1 M_sun, got {:.4} M_sun",
        mass_solar
    );

    println!("Total disk mass: {:.4} M_sun", mass_solar);
}

#[test]
fn viscous_timescale_reasonable() {
    let disk = GasDisk::mmsn();

    // Viscous timescale at 1 AU should be ~10^5 to 10^6 years for α = 10^(-3)
    let t_visc = disk.viscous_timescale(Length::from_au(1.0));
    let t_visc_yr = t_visc.to_years();

    assert!(
        t_visc_yr > 1e4 && t_visc_yr < 1e7,
        "t_visc at 1 AU: expected ~10^5-10^6 yr, got {:.2e} yr",
        t_visc_yr
    );

    println!("Viscous timescale at 1 AU: {:.2e} yr", t_visc_yr);
}

#[test]
fn mean_free_path_reasonable() {
    let disk = GasDisk::mmsn();

    // Mean free path at 1 AU should be ~1-100 cm
    let mfp = disk.mean_free_path(Length::from_au(1.0));

    assert!(
        mfp.to_cm() > 0.1 && mfp.to_cm() < 1000.0,
        "Mean free path at 1 AU: expected ~1-100 cm, got {:.1} cm",
        mfp.to_cm()
    );

    println!("Mean free path at 1 AU: {:.1} cm", mfp.to_cm());
}

#[test]
fn analytical_pressure_gradient_matches_expected() {
    let disk = GasDisk::mmsn();

    // For MMSN with p=1, q=0.5:
    // d ln P / d ln r = -(p + (3+q)/2) = -(1 + 1.75) = -2.75
    let d_ln_p = disk.pressure_gradient_log(Length::from_au(1.0));
    assert_relative_eq!(d_ln_p, -2.75, max_relative = 1e-10);

    // Should be constant across radius for power-law disk
    let d_ln_p_5au = disk.pressure_gradient_log(Length::from_au(5.0));
    assert_relative_eq!(d_ln_p, d_ln_p_5au, max_relative = 1e-10);
}

#[test]
fn tapered_disk_reduces_outer_density() {
    let disk_plain = GasDisk::mmsn();
    let disk_tapered = GasDisk::mmsn().with_standard_taper();

    // At 1 AU: taper should have minimal effect
    let sigma_1au_plain = disk_plain.surface_density(Length::from_au(1.0));
    let sigma_1au_tapered = disk_tapered.surface_density(Length::from_au(1.0));
    let ratio_1au = sigma_1au_tapered.to_grams_per_cm2() / sigma_1au_plain.to_grams_per_cm2();
    assert!(ratio_1au > 0.99, "At 1 AU, taper effect should be minimal");

    // At 40 AU (= r_c): taper reduces density to ~37%
    let sigma_40au_plain = disk_plain.surface_density(Length::from_au(40.0));
    let sigma_40au_tapered = disk_tapered.surface_density(Length::from_au(40.0));
    let ratio_40au = sigma_40au_tapered.to_grams_per_cm2() / sigma_40au_plain.to_grams_per_cm2();
    assert!(
        (ratio_40au - 0.368).abs() < 0.01,
        "At r_c, exp(-1) ≈ 0.368, got {}",
        ratio_40au
    );

    // At 100 AU (outer edge): taper reduces density dramatically
    let sigma_100au_plain = disk_plain.surface_density(Length::from_au(100.0));
    let sigma_100au_tapered = disk_tapered.surface_density(Length::from_au(100.0));
    let ratio_100au = sigma_100au_tapered.to_grams_per_cm2() / sigma_100au_plain.to_grams_per_cm2();
    assert!(
        ratio_100au < 0.01,
        "At outer edge, taper should reduce density to < 1%, got {}",
        ratio_100au
    );

    println!("Taper effect at 1 AU: {:.3}", ratio_1au);
    println!("Taper effect at 40 AU (r_c): {:.3}", ratio_40au);
    println!("Taper effect at 100 AU: {:.3}", ratio_100au);
}

#[test]
fn tapered_disk_with_custom_radius() {
    let disk_plain = GasDisk::mmsn();
    let custom_r_c = Length::from_au(20.0);
    let disk_tapered = GasDisk::mmsn().with_taper(custom_r_c);

    // At r_c = 20 AU: should have exp(-1) ≈ 0.368 of plain disk
    let sigma_plain = disk_plain.surface_density(custom_r_c);
    let sigma_tapered = disk_tapered.surface_density(custom_r_c);
    let ratio = sigma_tapered.to_grams_per_cm2() / sigma_plain.to_grams_per_cm2();

    assert!(
        (ratio - 0.368).abs() < 0.01,
        "At custom r_c = 20 AU, expected exp(-1) ≈ 0.368, got {}",
        ratio
    );
}

#[test]
fn tapered_disk_creates_different_pressure_profile() {
    let disk_plain = GasDisk::mmsn();
    let disk_tapered = GasDisk::mmsn().with_standard_taper();

    // In a plain power-law disk, pressure decreases monotonically
    // In a tapered disk, the exponential cutoff changes the pressure gradient

    // Sample at several radii
    let radii = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0];

    // Calculate pressure gradients (ratio of successive pressures)
    let mut plain_gradients = Vec::new();
    let mut tapered_gradients = Vec::new();

    for i in 0..radii.len() - 1 {
        let r1 = Length::from_au(radii[i]);
        let r2 = Length::from_au(radii[i + 1]);

        let p1_plain = disk_plain.surface_density(r1).to_grams_per_cm2()
            * disk_plain.temperature(r1).to_kelvin();
        let p2_plain = disk_plain.surface_density(r2).to_grams_per_cm2()
            * disk_plain.temperature(r2).to_kelvin();
        plain_gradients.push(p2_plain / p1_plain);

        let p1_tapered = disk_tapered.surface_density(r1).to_grams_per_cm2()
            * disk_tapered.temperature(r1).to_kelvin();
        let p2_tapered = disk_tapered.surface_density(r2).to_grams_per_cm2()
            * disk_tapered.temperature(r2).to_kelvin();
        tapered_gradients.push(p2_tapered / p1_tapered);
    }

    // In the outer disk (beyond r_c), tapered disk should have steeper pressure drop
    let outer_idx = tapered_gradients.len() - 1; // 60-80 AU range
    assert!(
        tapered_gradients[outer_idx] < plain_gradients[outer_idx],
        "Tapered disk should have steeper pressure gradient in outer regions"
    );

    println!(
        "Plain disk gradient (60-80 AU): {:.3}",
        plain_gradients[outer_idx]
    );
    println!(
        "Tapered disk gradient (60-80 AU): {:.3}",
        tapered_gradients[outer_idx]
    );
}

#[test]
fn no_taper_by_default() {
    let disk = GasDisk::mmsn();
    assert!(
        disk.taper_radius.is_none(),
        "Default disk should have no taper"
    );
}

#[test]
fn standard_taper_is_40_percent_of_outer_radius() {
    let disk = GasDisk::mmsn().with_standard_taper();
    let r_c = disk.taper_radius.expect("Should have taper radius");
    let expected = 0.4 * disk.outer_radius.to_au();

    assert_relative_eq!(r_c.to_au(), expected, max_relative = 1e-10);
}
