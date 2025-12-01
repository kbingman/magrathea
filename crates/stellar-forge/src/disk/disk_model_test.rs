use approx::assert_relative_eq;

use crate::disk::disk_model::DiskModel;
use units::{Length, Mass, SurfaceDensity, Temperature};

/// Minimal test disk implementing only required methods.
struct TestDisk {
    stellar_mass: Mass,
    sigma_0: f64,
    sigma_exponent: f64,
    temp_0: f64,
    temp_exponent: f64,
    alpha: f64,
    r_0: f64,
    inner_radius: Length,
    outer_radius: Length,
}

impl TestDisk {
    fn mmsn() -> Self {
        Self {
            stellar_mass: Mass::from_solar_masses(1.0),
            sigma_0: 1700.0, // g/cm² at 1 AU
            sigma_exponent: 1.0,
            temp_0: 280.0, // K at 1 AU
            temp_exponent: 0.5,
            alpha: 1e-3,
            r_0: 1.496e13, // 1 AU in cm
            inner_radius: Length::from_au(0.1),
            outer_radius: Length::from_au(100.0),
        }
    }
}

impl DiskModel for TestDisk {
    fn surface_density(&self, r: Length) -> SurfaceDensity {
        let ratio = r.to_cm() / self.r_0;
        SurfaceDensity::from_grams_per_cm2(self.sigma_0 * ratio.powf(-self.sigma_exponent))
    }

    fn temperature(&self, r: Length) -> Temperature {
        let ratio = r.to_cm() / self.r_0;
        Temperature::from_kelvin(self.temp_0 * ratio.powf(-self.temp_exponent))
    }

    fn stellar_mass(&self) -> Mass {
        self.stellar_mass
    }

    fn alpha(&self) -> f64 {
        self.alpha
    }

    fn inner_radius(&self) -> Length {
        self.inner_radius
    }

    fn outer_radius(&self) -> Length {
        self.outer_radius
    }
}

#[test]
fn surface_density_at_1au() {
    let disk = TestDisk::mmsn();
    let sigma = disk.surface_density(Length::from_au(1.0));
    assert_relative_eq!(sigma.to_grams_per_cm2(), 1700.0, max_relative = 1e-10);
}

#[test]
fn temperature_at_1au() {
    let disk = TestDisk::mmsn();
    let t = disk.temperature(Length::from_au(1.0));
    assert_relative_eq!(t.to_kelvin(), 280.0, max_relative = 1e-10);
}

#[test]
fn orbital_period_at_1au() {
    let disk = TestDisk::mmsn();
    let period = disk.orbital_period(Length::from_au(1.0));
    // Should be ~1 year
    assert_relative_eq!(period.to_years(), 1.0, max_relative = 0.01);
}

#[test]
fn aspect_ratio_reasonable() {
    let disk = TestDisk::mmsn();
    let h_r = disk.aspect_ratio(Length::from_au(1.0));
    // Should be ~0.03-0.05 at 1 AU
    assert!(h_r > 0.02 && h_r < 0.07, "h/r at 1 AU: {:.4}", h_r);
}

#[test]
fn pressure_gradient_parameter_reasonable() {
    let disk = TestDisk::mmsn();
    let eta = disk.pressure_gradient_parameter(Length::from_au(1.0));
    // Should be ~0.002-0.005
    assert!(eta > 0.001 && eta < 0.01, "η at 1 AU: {:.5}", eta);
}

#[test]
fn numerical_pressure_gradient_reasonable() {
    let disk = TestDisk::mmsn();
    // For MMSN: d ln P / d ln r = -(p + (3+q)/2) = -(1 + 1.75) = -2.75
    let d_ln_p = disk.pressure_gradient_log(Length::from_au(1.0));
    assert_relative_eq!(d_ln_p, -2.75, max_relative = 0.01);
}

#[test]
fn viscous_timescale_reasonable() {
    let disk = TestDisk::mmsn();
    let t_visc = disk.viscous_timescale(Length::from_au(1.0));
    // Should be ~10^5 to 10^6 years for α = 10^(-3)
    let t_yr = t_visc.to_years();
    assert!(t_yr > 1e4 && t_yr < 1e7, "t_visc at 1 AU: {:.2e} yr", t_yr);
}

#[test]
fn is_valid_radius_works() {
    let disk = TestDisk::mmsn();
    assert!(!disk.is_valid_radius(Length::from_au(0.05)));
    assert!(disk.is_valid_radius(Length::from_au(1.0)));
    assert!(disk.is_valid_radius(Length::from_au(50.0)));
    assert!(!disk.is_valid_radius(Length::from_au(150.0)));
}
