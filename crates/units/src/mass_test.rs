mod tests {
    use approx::assert_relative_eq;

    use crate::mass::{Mass, SOLAR_MASS_G};

    #[test]
    fn test_mass_conversions() {
        // Test solar masses to grams
        let mass_sm = Mass::from_solar_masses(1.0);
        assert_relative_eq!(mass_sm.to_grams(), SOLAR_MASS_G);

        // Test grams to solar masses
        let mass_g = Mass::from_grams(SOLAR_MASS_G);
        assert_relative_eq!(mass_g.to_solar_masses(), 1.0);

        // Test round trip
        let original = 0.05; // A small star
        let mass = Mass::from_solar_masses(original);
        let g_value = mass.to_grams();
        let round_trip = Mass::from_grams(g_value).to_solar_masses();
        assert_relative_eq!(round_trip, original);
    }

    #[test]
    fn test_mass_arithmetic_operations() {
        let mass1 = Mass::from_solar_masses(2.0);
        let mass2 = Mass::from_solar_masses(1.5);

        // Test addition and subtraction
        assert_relative_eq!((mass1 + mass2).to_solar_masses(), 3.5);
        assert_relative_eq!((mass1 - mass2).to_solar_masses(), 0.5);

        // Test multiplication with f64
        let scaled = mass1 * 3.0;
        assert_relative_eq!(scaled.to_solar_masses(), 6.0);

        // Test division with f64
        let divided = mass1 / 4.0;
        assert_relative_eq!(divided.to_solar_masses(), 0.5);

        // Test commutative multiplication
        let mass = Mass::from_earth_masses(100.0);
        let commutative = 2.5 * mass;
        assert_relative_eq!(commutative.to_earth_masses(), 250.0);
    }
}
