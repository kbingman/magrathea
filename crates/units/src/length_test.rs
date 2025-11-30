mod tests {
    use approx::assert_relative_eq;

    use crate::length::{Length, AU_TO_CM};

    #[test]
    fn test_length_conversions() {
        // Test AU to cm conversion
        let length_au = Length::from_au(1.0);
        assert_relative_eq!(length_au.to_cm(), AU_TO_CM);

        // Test cm to AU conversion
        let length_cm = Length::from_cm(AU_TO_CM);
        assert_relative_eq!(length_cm.to_au(), 1.0);

        // Test round trip
        let original = 5.7;
        let length = Length::from_au(original);
        let cm_value = length.to_cm();
        let round_trip = Length::from_cm(cm_value).to_au();
        assert_relative_eq!(round_trip, original);
    }

    #[test]
    fn test_length_arithmetic_operations() {
        let length1 = Length::from_au(5.0);
        let length2 = Length::from_au(3.0);

        // Test addition and subtraction
        assert_relative_eq!((length1 + length2).to_au(), 8.0);
        assert_relative_eq!((length1 - length2).to_au(), 2.0);

        // Test multiplication with f64
        let scaled = length1 * 2.0;
        assert_relative_eq!(scaled.to_au(), 10.0);

        // Test division with f64
        let divided = length1 / 2.0;
        assert_relative_eq!(divided.to_au(), 2.5);

        // Test commutative multiplication
        let commutative = 1.5 * length1;
        assert_relative_eq!(commutative.to_au(), 7.5);
    }

    #[test]
    fn test_length_min_max() {
        let length1 = Length::from_au(5.0);
        let length2 = Length::from_au(3.0);
        let length3 = Length::from_au(7.0);

        // Test min
        assert_relative_eq!(length1.min(length2).to_au(), 3.0);
        assert_relative_eq!(length2.min(length1).to_au(), 3.0);
        assert_relative_eq!(length1.min(length3).to_au(), 5.0);

        // Test max
        assert_relative_eq!(length1.max(length2).to_au(), 5.0);
        assert_relative_eq!(length2.max(length1).to_au(), 5.0);
        assert_relative_eq!(length1.max(length3).to_au(), 7.0);
    }
}
