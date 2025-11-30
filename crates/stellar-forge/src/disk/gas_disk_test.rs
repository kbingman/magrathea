use crate::disk::gas_disk::GasDisk;
use units::Length;

/// Helper to check if two f64 values are approximately equal.
fn approx_eq(a: f64, b: f64, rel_tol: f64) -> bool {
    let diff = (a - b).abs();
    let max = a.abs().max(b.abs());
    if max == 0.0 {
        diff < 1e-15
    } else {
        diff / max < rel_tol
    }
}

#[test]
fn mmsn_surface_density_at_1au() {
    let disk = GasDisk::mmsn();
    let sigma = disk.surface_density(Length::from_au(1.0));

    // Should be 1700 g/cm² at 1 AU by construction
    assert!(
        approx_eq(sigma.to_grams_per_cm2(), 1700.0, 1e-10),
        "Σ at 1 AU: expected 1700, got {}",
        sigma.to_grams_per_cm2()
    );
}

#[test]
fn mmsn_surface_density_at_5au() {
    let disk = GasDisk::mmsn();
    let sigma = disk.surface_density(Length::from_au(5.0));

    // Σ(5 AU) = 1700 × 5^(-1) = 340 g/cm²
    let expected = 1700.0 / 5.0;
    assert!(
        approx_eq(sigma.to_grams_per_cm2(), expected, 1e-10),
        "Σ at 5 AU: expected {}, got {}",
        expected,
        sigma.to_grams_per_cm2()
    );
}

#[test]
fn mmsn_temperature_at_1au() {
    let disk = GasDisk::mmsn();
    let t = disk.temperature(Length::from_au(1.0));

    // Should be 280 K at 1 AU by construction
    assert!(
        approx_eq(t.to_kelvin(), 280.0, 1e-10),
        "T at 1 AU: expected 280, got {}",
        t.to_kelvin()
    );
}

#[test]
fn mmsn_temperature_at_5au() {
    let disk = GasDisk::mmsn();
    let t = disk.temperature(Length::from_au(5.0));

    // T(5 AU) = 280 × 5^(-0.5) ≈ 125 K
    let expected = 280.0 / 5.0_f64.sqrt();
    assert!(
        approx_eq(t.to_kelvin(), expected, 1e-10),
        "T at 5 AU: expected {:.1}, got {:.1}",
        expected,
        t.to_kelvin()
    );
}

#[test]
fn orbital_period_at_1au() {
    let disk = GasDisk::mmsn();
    let period = disk.orbital_period(Length::from_au(1.0));

    // Should be ~1 year
    assert!(
        approx_eq(period.to_years(), 1.0, 0.01),
        "Period at 1 AU: expected ~1 yr, got {:.3} yr",
        period.to_years()
    );
}

#[test]
fn orbital_period_at_5au() {
    let disk = GasDisk::mmsn();
    let period = disk.orbital_period(Length::from_au(5.0));

    // Kepler's third law: P ∝ a^(3/2), so P(5 AU) ≈ 11.2 yr
    let expected = 5.0_f64.powf(1.5);
    assert!(
        approx_eq(period.to_years(), expected, 0.01),
        "Period at 5 AU: expected {:.1} yr, got {:.1} yr",
        expected,
        period.to_years()
    );
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
