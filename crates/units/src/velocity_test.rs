mod tests {
    use approx::assert_relative_eq;

    use crate::velocity::{
        circular_orbital_velocity, Velocity, AU_YEAR_TO_CM_SEC, EARTH_ORBITAL_VELOCITY_AU_PER_YEAR,
    };

    #[test]
    fn test_velocity_conversions() {
        // Test AU/year to cm/sec conversion
        let vel_au_year = Velocity::from_au_per_year(1.0);
        assert_relative_eq!(vel_au_year.to_cm_per_sec(), AU_YEAR_TO_CM_SEC);

        // Test cm/sec to AU/year conversion
        let vel_cm_sec = Velocity::from_cm_per_sec(AU_YEAR_TO_CM_SEC);
        assert_relative_eq!(vel_cm_sec.to_au_per_year(), 1.0);

        // Test round trip
        let original = 0.75; // Some orbital velocity
        let velocity = Velocity::from_au_per_year(original);
        let cm_sec_value = velocity.to_cm_per_sec();
        let round_trip = Velocity::from_cm_per_sec(cm_sec_value).to_au_per_year();
        assert_relative_eq!(round_trip, original);
    }

    #[test]
    fn test_velocity_operations() {
        // Test addition
        let a = Velocity::from_au_per_year(1.0);
        let b = Velocity::from_au_per_year(2.0);
        let sum = a + b;
        assert_relative_eq!(sum.to_au_per_year(), 3.0);

        // Test subtraction
        let diff = b - a;
        assert_relative_eq!(diff.to_au_per_year(), 1.0);

        // Test multiplication
        let mult = a * 3.0;
        assert_relative_eq!(mult.to_au_per_year(), 3.0);

        // Test division
        let div = b / 2.0;
        assert_relative_eq!(div.to_au_per_year(), 1.0);
    }

    #[test]
    fn test_circular_orbital_velocity() {
        // Test Earth's orbital velocity around 1 solar mass at 1 AU
        let earth_velocity = circular_orbital_velocity(1.0, 1.0);
        assert_relative_eq!(earth_velocity, EARTH_ORBITAL_VELOCITY_AU_PER_YEAR);

        // Test Jupiter's orbital velocity around 1 solar mass at 5.2 AU
        let jupiter_velocity = circular_orbital_velocity(1.0, 5.2);
        let expected_jupiter_velocity = EARTH_ORBITAL_VELOCITY_AU_PER_YEAR / (5.2_f32).sqrt();
        assert_relative_eq!(jupiter_velocity, expected_jupiter_velocity, epsilon = 0.1);

        // Test Mercury's orbital velocity around 1 solar mass at 0.39 AU
        let mercury_velocity = circular_orbital_velocity(1.0, 0.39);
        let expected_mercury_velocity = EARTH_ORBITAL_VELOCITY_AU_PER_YEAR / (0.39_f32).sqrt();
        assert_relative_eq!(mercury_velocity, expected_mercury_velocity, epsilon = 0.1);

        // Test scaling with stellar mass - 2 solar mass star
        let massive_star_velocity = circular_orbital_velocity(2.0, 1.0);
        let expected_massive_velocity = EARTH_ORBITAL_VELOCITY_AU_PER_YEAR * (2.0_f32).sqrt();
        assert_relative_eq!(massive_star_velocity, expected_massive_velocity);
    }

    #[test]
    fn test_circular_orbital_velocity_real_world_values() {
        // Test against known real-world orbital velocities in AU/year
        // Earth: ~29.78 km/s = ~6.28 AU/year

        // Earth: should give ~6.28 AU/year (TAU ≈ 2π)
        let earth_v = circular_orbital_velocity(1.0, 1.0);
        assert_relative_eq!(earth_v, std::f64::consts::TAU as f32, epsilon = 0.01);

        // Jupiter: should give ~2.75 AU/year (6.28 / sqrt(5.2) ≈ 2.75)
        let jupiter_v = circular_orbital_velocity(1.0, 5.2);
        assert_relative_eq!(jupiter_v, 2.75, epsilon = 0.1);

        // Mercury: should give ~10.07 AU/year (6.28 / sqrt(0.39) ≈ 10.07)
        let mercury_v = circular_orbital_velocity(1.0, 0.39);
        assert_relative_eq!(mercury_v, 10.07, epsilon = 0.1);
    }
}
