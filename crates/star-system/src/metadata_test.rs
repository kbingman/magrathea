use crate::{GenerationMethod, SystemArchitecture, SystemMetadata};

#[test]
fn test_catalog_name_format() {
    let meta = SystemMetadata::new_random(GenerationMethod::Statistical, SystemArchitecture::Mixed);
    let name = meta.catalog_name();

    // Format: XX-0000
    assert_eq!(name.len(), 7);
    assert!(name.chars().nth(0).unwrap().is_ascii_uppercase());
    assert!(name.chars().nth(1).unwrap().is_ascii_uppercase());
    assert_eq!(name.chars().nth(2).unwrap(), '-');
    assert!(name[3..].chars().all(|c| c.is_ascii_digit()));
}

#[test]
fn test_catalog_name_deterministic() {
    let meta1 = SystemMetadata::from_seed_name(
        "test",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );
    let meta2 = SystemMetadata::from_seed_name(
        "test",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );

    assert_eq!(meta1.catalog_name(), meta2.catalog_name());
}

#[test]
fn test_seed_deterministic() {
    let meta1 = SystemMetadata::from_seed_name(
        "test",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );
    let meta2 = SystemMetadata::from_seed_name(
        "test",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );

    assert_eq!(meta1.seed(), meta2.seed());
}

#[test]
fn test_display_name_without_proper_name() {
    let meta = SystemMetadata::new_random(GenerationMethod::Statistical, SystemArchitecture::Mixed);
    assert_eq!(meta.display_name(), meta.catalog_name());
}

#[test]
fn test_display_name_with_proper_name() {
    let meta = SystemMetadata::new_random(GenerationMethod::Statistical, SystemArchitecture::Mixed)
        .with_name("Cygnus Prime");
    assert_eq!(meta.display_name(), "Cygnus Prime");
    // Catalog name still available
    assert_ne!(meta.catalog_name(), "Cygnus Prime");
}

#[test]
fn test_different_seeds_different_names() {
    let meta1 = SystemMetadata::from_seed_name(
        "alpha",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );
    let meta2 = SystemMetadata::from_seed_name(
        "beta",
        GenerationMethod::Statistical,
        SystemArchitecture::Mixed,
    );

    assert_ne!(meta1.catalog_name(), meta2.catalog_name());
}
