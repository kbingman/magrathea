mod tests {
    use approx::assert_relative_eq;

    use crate::temperature::Temperature;

    #[test]
    fn test_temperature_conversions() {
        // Test Kelvin to Celsius
        let freezing = Temperature::from_kelvin(273.15);
        assert_relative_eq!(freezing.to_celsius(), 0.0, epsilon = 0.01);

        // Test Celsius to Kelvin round trip
        let temp_c = Temperature::from_celsius(100.0);
        assert_relative_eq!(temp_c.to_kelvin(), 373.15, epsilon = 0.01);

        // Test Fahrenheit conversions
        let temp_f = Temperature::from_fahrenheit(32.0);
        assert_relative_eq!(temp_f.to_kelvin(), 273.15, epsilon = 0.01);
        assert_relative_eq!(temp_f.to_celsius(), 0.0, epsilon = 0.01);

        // 68°F should be about 20°C
        let room_f = Temperature::from_fahrenheit(68.0);
        assert_relative_eq!(room_f.to_celsius(), 20.0, epsilon = 0.5);
    }

    #[test]
    fn test_temperature_round_trips() {
        // Kelvin -> Celsius -> Kelvin
        let original_k = 300.0;
        let temp = Temperature::from_kelvin(original_k);
        let celsius = temp.to_celsius();
        let round_trip = Temperature::from_celsius(celsius);
        assert_relative_eq!(round_trip.to_kelvin(), original_k, epsilon = 0.01);

        // Kelvin -> Fahrenheit -> Kelvin
        let temp2 = Temperature::from_kelvin(273.15);
        let fahrenheit = temp2.to_fahrenheit();
        let round_trip2 = Temperature::from_fahrenheit(fahrenheit);
        assert_relative_eq!(round_trip2.to_kelvin(), 273.15, epsilon = 0.01);
    }

    #[test]
    fn test_temperature_arithmetic() {
        let temp1 = Temperature::from_kelvin(300.0);
        let temp2 = Temperature::from_kelvin(50.0);

        // Test addition and subtraction
        assert_relative_eq!((temp1 + temp2).to_kelvin(), 350.0);
        assert_relative_eq!((temp1 - temp2).to_kelvin(), 250.0);

        // Test multiplication and division
        let doubled = temp1 * 2.0;
        assert_relative_eq!(doubled.to_kelvin(), 600.0);

        let halved = temp1 / 2.0;
        assert_relative_eq!(halved.to_kelvin(), 150.0);
    }

    #[test]
    fn test_temperature_constants() {
        let freezing = Temperature::water_freezing();
        assert_relative_eq!(freezing.to_kelvin(), 273.15, epsilon = 0.01);
        assert_relative_eq!(freezing.to_celsius(), 0.0, epsilon = 0.01);

        let boiling = Temperature::water_boiling();
        assert_relative_eq!(boiling.to_kelvin(), 373.15, epsilon = 0.01);
        assert_relative_eq!(boiling.to_celsius(), 100.0, epsilon = 0.01);
    }

    #[test]
    fn test_unit_conversion_consistency() {
        // 0°C = 32°F = 273.15 K
        let freezing_c = Temperature::from_celsius(0.0);
        let freezing_f = Temperature::from_fahrenheit(32.0);
        assert_relative_eq!(
            freezing_c.to_kelvin(),
            freezing_f.to_kelvin(),
            epsilon = 0.01
        );

        // 100°C = 212°F = 373.15 K
        let boiling_c = Temperature::from_celsius(100.0);
        let boiling_f = Temperature::from_fahrenheit(212.0);
        assert_relative_eq!(boiling_c.to_kelvin(), boiling_f.to_kelvin(), epsilon = 0.01);
    }
}
