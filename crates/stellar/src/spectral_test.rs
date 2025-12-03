use crate::spectral::{LuminosityClass, SpectralType, VariabilityType};

// ============================================================================
// SpectralType Display tests
// ============================================================================

#[test]
fn spectral_type_display_main_sequence() {
    assert_eq!(format!("{}", SpectralType::O), "O");
    assert_eq!(format!("{}", SpectralType::B), "B");
    assert_eq!(format!("{}", SpectralType::A), "A");
    assert_eq!(format!("{}", SpectralType::F), "F");
    assert_eq!(format!("{}", SpectralType::G), "G");
    assert_eq!(format!("{}", SpectralType::K), "K");
    assert_eq!(format!("{}", SpectralType::M), "M");
}

#[test]
fn spectral_type_display_brown_dwarfs() {
    assert_eq!(format!("{}", SpectralType::L), "L");
    assert_eq!(format!("{}", SpectralType::T), "T");
    assert_eq!(format!("{}", SpectralType::Y), "Y");
}

#[test]
fn spectral_type_display_remnants() {
    assert_eq!(format!("{}", SpectralType::D), "D");
    assert_eq!(format!("{}", SpectralType::N), "N");
    assert_eq!(format!("{}", SpectralType::BH), "BH");
}

#[test]
fn spectral_type_display_rare_types() {
    assert_eq!(format!("{}", SpectralType::S), "S");
    assert_eq!(format!("{}", SpectralType::C), "C");
    assert_eq!(format!("{}", SpectralType::R), "R");
    assert_eq!(format!("{}", SpectralType::Q), "Q");
}

// ============================================================================
// LuminosityClass Display tests
// ============================================================================

#[test]
fn luminosity_class_display_supergiants() {
    assert_eq!(format!("{}", LuminosityClass::IAPLUS), "Ia+");
    assert_eq!(format!("{}", LuminosityClass::IA), "Ia");
    assert_eq!(format!("{}", LuminosityClass::IB), "Ib");
}

#[test]
fn luminosity_class_display_giants() {
    assert_eq!(format!("{}", LuminosityClass::II), "II");
    assert_eq!(format!("{}", LuminosityClass::III), "III");
    assert_eq!(format!("{}", LuminosityClass::IV), "IV");
}

#[test]
fn luminosity_class_display_dwarfs() {
    assert_eq!(format!("{}", LuminosityClass::V), "V");
    assert_eq!(format!("{}", LuminosityClass::VI), "VI");
    assert_eq!(format!("{}", LuminosityClass::D), "D");
    assert_eq!(format!("{}", LuminosityClass::BD), "BD");
    assert_eq!(format!("{}", LuminosityClass::NS), "NS");
}

// ============================================================================
// VariabilityType::determine_variability tests
// ============================================================================

#[test]
fn variability_cepheid_supergiants_in_instability_strip() {
    // Cepheid: supergiant (Ia/Ib), mass > 5, temp 5000-7000K
    let variability = VariabilityType::determine_variability(10.0, 6000.0, LuminosityClass::IA);
    assert_eq!(variability, VariabilityType::Cepheid);

    let variability = VariabilityType::determine_variability(8.0, 5500.0, LuminosityClass::IB);
    assert_eq!(variability, VariabilityType::Cepheid);
}

#[test]
fn variability_rr_lyrae_horizontal_branch() {
    // RR Lyrae: temp 6000-7500K (any luminosity class not matching Cepheid)
    let variability = VariabilityType::determine_variability(1.0, 6500.0, LuminosityClass::V);
    assert_eq!(variability, VariabilityType::RRLyrae);

    let variability = VariabilityType::determine_variability(2.0, 7000.0, LuminosityClass::IV);
    assert_eq!(variability, VariabilityType::RRLyrae);
}

#[test]
fn variability_mira_cool_giants() {
    // Mira: cool giants (III), temp < 3500K
    let variability = VariabilityType::determine_variability(2.0, 3000.0, LuminosityClass::III);
    assert_eq!(variability, VariabilityType::Mira);
}

#[test]
fn variability_beta_cephei_hot_massive() {
    // Beta Cephei: mass > 8, temp > 20000K
    let variability = VariabilityType::determine_variability(12.0, 25000.0, LuminosityClass::V);
    assert_eq!(variability, VariabilityType::Beta);
}

#[test]
fn variability_delta_scuti_af_dwarfs() {
    // Delta Scuti: main sequence (V), temp 6500-8500K
    // Note: RR Lyrae check (6000-7500K) takes precedence, so test at high end
    let variability = VariabilityType::determine_variability(1.8, 8000.0, LuminosityClass::V);
    assert_eq!(variability, VariabilityType::Delta);
}

#[test]
fn variability_flare_m_dwarfs() {
    // Flare stars: main sequence (V), temp < 3500K
    let variability = VariabilityType::determine_variability(0.3, 3000.0, LuminosityClass::V);
    assert_eq!(variability, VariabilityType::Flare);
}

#[test]
fn variability_none_default() {
    // Non-variable: anything not matching other patterns
    let variability = VariabilityType::determine_variability(1.0, 5800.0, LuminosityClass::V);
    assert_eq!(variability, VariabilityType::None);

    // Giant but not cool enough for Mira
    let variability = VariabilityType::determine_variability(2.0, 4500.0, LuminosityClass::III);
    assert_eq!(variability, VariabilityType::None);
}
