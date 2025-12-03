# Occurrence Rate Calibration for Statistical Planet Generator

## Summary of Literature Values

Based on extensive literature review (Fressin et al. 2013, Petigura et al. 2013, Howard et al. 2012, Burke et al. 2015, Dressing & Charbonneau 2015, Mulders et al. 2015, 2018, Hsu et al. 2019, Fulton et al. 2017, and others), here are the calibrated occurrence rates.

---

## 1. Planet Occurrence by Size (FGK Stars, P < 100 days)

These are planets per star in each radius bin:

| Radius Range | Name | Occurrence Rate | Source |
|--------------|------|-----------------|--------|
| 0.5–1.0 R⊕ | Sub-Earth | ~5–10% | Fressin+ 2013, extrapolation |
| 1.0–1.25 R⊕ | Earth-size | ~15% (14.9 ± 2.4% at P < 50d) | Fressin+ 2013 |
| 1.25–2.0 R⊕ | Super-Earth | ~20–25% | Fressin+ 2013, Petigura+ 2013 |
| 2.0–4.0 R⊕ | Sub-Neptune | ~25–30% | Most common class |
| 4.0–6.0 R⊕ | Neptune | ~5–8% | |
| 6.0–10 R⊕ | Sub-Saturn | ~2–3% | Rare, "Neptune desert" |
| 10–22 R⊕ | Jupiter | ~1–2% | |
| Hot Jupiters | P < 10d, R > 10 R⊕ | ~0.5–1% | Wright+ 2012 |

**Key finding**: Sub-Neptunes (2–4 R⊕) are the most common planet type. The distribution peaks around 2.4 R⊕ with a valley at 1.5–2.0 R⊕ (Fulton gap).

---

## 2. Planet Occurrence by Period (All Sizes Combined)

| Period Range | Planets per Star | Notes |
|--------------|------------------|-------|
| < 3 days | ~2% | Ultra-short period |
| 3–10 days | ~15% | Hot zone |
| 10–50 days | ~35% | Peak occurrence |
| 50–100 days | ~25% | Warm zone |
| 100–200 days | ~15% | |
| 200–400 days | ~10% | Near HZ for G stars |

**Key finding**: Occurrence rate increases with orbital period up to ~10 days, then flattens. This ~10 day break corresponds to the inner edge of planetary systems (Mulders+ 2018).

---

## 3. Spectral Type Dependence

Sub-Neptunes (R < 4 R⊕) are **~3× more common** around M dwarfs than FGK stars (Mulders+ 2015, 2018).

| Stellar Type | Sub-Neptune Occurrence | Giant Planet Occurrence |
|--------------|------------------------|------------------------|
| M dwarfs | ~2.5 planets/star | ~3% |
| K dwarfs | ~1.5 planets/star | ~5% |
| G dwarfs | ~1.0 planets/star | ~10% |
| F dwarfs | ~0.7 planets/star | ~10–15% |
| A dwarfs | ~0.3 planets/star | Poorly constrained |

**Key insight**: M dwarfs prefer compact, multi-planet systems with small planets. FGK stars have more diverse architectures including giants.

---

## 4. Metallicity Dependence

### Giants (R > 6 R⊕)
**Strong correlation**: Fischer & Valenti 2005 established:
```
P(giant) ∝ 10^(2 × [Fe/H])
```

| [Fe/H] | Giant Planet Rate |
|--------|-------------------|
| -0.5 | ~1% |
| 0.0 | ~5% |
| +0.3 | ~15–20% |
| +0.5 | ~25–30% |

### Small Planets (R < 4 R⊕)
**Weak or no correlation**: Small planets occur around metal-poor and metal-rich stars at similar rates.

| [Fe/H] | Super-Earth/Sub-Neptune Rate |
|--------|------------------------------|
| -0.5 to -0.3 | ~20–25% |
| -0.3 to 0.0 | ~25–30% |
| 0.0 to +0.3 | ~30–35% |
| +0.3 to +0.5 | ~30–35% |

---

## 5. Integrated Statistics

From Mulders+ 2018 and Hsu+ 2019:

- **At least 42% of Sun-like stars** have coplanar multi-planet systems with 7+ planets
- **~100% of stars** may have at least one planet within 1 AU (model-dependent)
- **η⊕ (Earth-like in HZ)**: 36 ± 14% for Sun-like stars (Mulders+ 2018); upper limit ~0.27 for FGK (Hsu+ 2019)
- **~30% of Sun-like stars** have Kepler-like planetary systems (Zhu+ 2018)

---

## 6. Outer System Occurrence Rates (Beyond Kepler)

Kepler's transit method is biased toward short periods (P < 400 days, a < 1.2 AU). For outer system planets, we rely on radial velocity surveys, direct imaging, and microlensing.

### 6.1 Radial Velocity Constraints (1–10 AU)

From the HARPS/CORALIE survey (Mayor+ 2011) and Cumming+ 2008:

| Planet Type | Semi-major Axis | Occurrence (FGK) | Notes |
|-------------|-----------------|------------------|-------|
| Cold Jupiters | 1–5 AU | ~10–15% | Peak at 2–3 AU |
| Cold Jupiters | 5–10 AU | ~5–8% | Declining with distance |
| Saturn-mass | 1–10 AU | ~5–10% | Less constrained |
| Ice Giants | 1–10 AU | ~10–20% | Detection limited |

**Key findings:**
- Giant planet occurrence rises steeply from ~3% at 0.1 AU to ~15% at 2–3 AU
- Beyond the snow line peak, occurrence declines as ~a^(-0.6) (Cumming+ 2008)
- ~14% of FGK stars have a giant planet (M > 50 M⊕) within 20 AU

### 6.2 Cold Giant Distribution (Fernandes+ 2019)

Meta-analysis of RV surveys gives the semi-major axis distribution:

```
dN/d(log a) ∝ a^β  where β ≈ 0.4 for a < 3 AU, β ≈ -0.5 for a > 3 AU
```

| Distance | Relative Occurrence | Cumulative |
|----------|---------------------|------------|
| 0.5–1 AU | 0.15 | 15% |
| 1–2 AU | 0.25 | 40% |
| 2–3 AU | 0.30 | 70% (peak) |
| 3–5 AU | 0.18 | 88% |
| 5–10 AU | 0.10 | 98% |
| >10 AU | 0.02 | 100% |

### 6.3 Microlensing Constraints (1–10 AU)

Microlensing is sensitive to planets near the Einstein ring (~2–4 AU for typical lens distances).

From Cassan+ 2012 and Suzuki+ 2016:

| Planet Mass | Occurrence per Star | Notes |
|-------------|---------------------|-------|
| Jupiter-mass (0.3–10 M_J) | 17 ± 6% | Per dex of mass |
| Neptune-mass (10–30 M⊕) | 52 ± 22% | Very common |
| Super-Earths (5–10 M⊕) | ~60% | Extrapolated |

**Key insight**: Microlensing suggests **cold Neptunes are very common** (~50% of stars), more so than Kepler's inner system rates would suggest. This supports abundant ice giant formation beyond the snow line.

### 6.4 Direct Imaging Constraints (10–1000 AU)

Direct imaging probes wide separations but only massive planets (>1 M_J). From GPIES (Nielsen+ 2019) and SPHERE/SHINE (Vigan+ 2021):

| Separation | Planet Mass | Occurrence | Notes |
|------------|-------------|------------|-------|
| 10–100 AU | 5–13 M_J | 0.8 ± 0.5% | GPIES |
| 10–100 AU | 1–20 M_J | 1–3% | Combined surveys |
| 100–300 AU | 5–13 M_J | <1% | Very rare |
| 10–1000 AU | >5 M_J | 1–5% | Total wide companions |

**Implications for generation:**
- Wide (50–200 AU) massive companions exist but are rare (~1–3%)
- These may form via disk instability or capture
- Include as rare "wide companion" population

### 6.5 Ice Giant Occurrence

Combining RV, microlensing, and Solar System analogy:

| Zone | Distance (Sun-like) | Ice Giant Occurrence | Notes |
|------|---------------------|----------------------|-------|
| Inner ice zone | 5–15 AU | ~15–25% | Like failed Saturns |
| Outer ice zone | 15–40 AU | ~10–20% | Uranus/Neptune analogs |
| Combined | 5–40 AU | ~25–40% | At least one ice giant |

**Solar System context**: Our system has 2 ice giants at 19 and 30 AU. Microlensing suggests this is common, not exceptional.

### 6.6 Recommended Outer System Zones

For a Sun-like star (L = 1 L☉, snow line ≈ 2.7 AU):

| Zone | Distance | Scaled by Snow Line | Primary Population |
|------|----------|---------------------|-------------------|
| Jupiter zone | 3–10 AU | 1–4× SL | Gas giants |
| Saturn zone | 8–18 AU | 3–7× SL | Gas giants, ice giants |
| Ice giant zone | 15–40 AU | 6–15× SL | Neptunes, Uranuses |
| Kuiper analog | 30–100 AU | 11–37× SL | Dwarf planets, debris |
| Wide companions | 50–500 AU | 18–185× SL | Rare massive planets |

### 6.7 Outer System Implementation

```rust
/// Outer system occurrence rates by zone
/// Probabilities are for FGK stars at solar metallicity
struct OuterSystemRates {
    // Cold giants: 1.5–6× snow line (~4–16 AU for Sun)
    // Base rate ~15%, boosted by metallicity
    cold_giant_base: f64,       // 0.15
    cold_giant_boost: f64,      // 10^(2×[Fe/H])

    // Ice giants: 5–12× snow line (~14–32 AU for Sun)
    // ~30% of systems, weak metallicity dependence
    ice_giant_rate: f64,        // 0.30

    // Wide companions: 50–300 AU
    // ~2% of systems, massive planets only
    wide_companion_rate: f64,   // 0.02
    wide_companion_min_mass: f64, // 100.0 (Earth masses, ~0.3 M_J)
}

/// Number of ice giants follows: 40% have 1, 35% have 2, 20% have 0, 5% have 3
fn sample_ice_giant_count(rng: &mut ChaChaRng) -> usize {
    match rng.random::<f64>() {
        x if x < 0.20 => 0,
        x if x < 0.60 => 1,
        x if x < 0.95 => 2,
        _ => 3,
    }
}
```

---

## 7. Architecture Statistics

From Zhu & Dong 2021 review:

| Architecture | Fraction of Stars | Characteristics |
|--------------|-------------------|-----------------|
| Compact Multi | ~30–40% | 4–8 planets, P < 100d, mutual i < 3° |
| Mixed | ~20–30% | Inner small + outer giant |
| Giant-dominated | ~5–10% | Hot/cold Jupiters, few companions |
| Sparse/Empty | ~30–40% | 0–1 detected planets |

### Hot Jupiter Systems
- Hot Jupiters are "lonely": systems with HJ have fewer companions
- ~0.5–1% of Sun-like stars have a hot Jupiter

### Period Ratio Distribution
- Excess near (but not at) mean-motion resonances
- Peak at P2/P1 ≈ 1.5–2.0

---

## 8. Recommended Implementation

### Mass Sampling Weights
Convert radius rates to mass using M-R relations. For FGK stars, [Fe/H] = 0:

```rust
/// Occurrence weights by mass bin (planets per star in each bin)
/// Normalized to sum to ~1.0 for relative sampling
const MASS_BIN_WEIGHTS: [(f64, f64, f64); 7] = [
    // (min_mass_earth, max_mass_earth, relative_weight)
    (0.01, 0.5,  0.08),  // Sub-Earth: less common, detection limited
    (0.5,  2.0,  0.22),  // Earth-sized: common
    (2.0,  5.0,  0.30),  // Super-Earth/Transitional: most common
    (5.0,  20.0, 0.25),  // Sub-Neptune/Volatile: very common
    (20.0, 50.0, 0.08),  // Neptune-class: less common
    (50.0, 160.0, 0.04), // Sub-Saturn: rare (Neptune desert)
    (160.0, 1000.0, 0.03), // Giant: rare but metallicity-boosted
];

/// Metallicity boost for giants
fn giant_metallicity_factor(feh: f64) -> f64 {
    10.0_f64.powf(2.0 * feh).clamp(0.1, 10.0)
}
```

### Period Sampling by Class
```rust
/// Period distribution parameters (days) by planet class
/// (inner_edge, peak, outer_characteristic, outer_falloff_power)
const PERIOD_PARAMS: [(PlanetClass, f64, f64, f64, f64); 4] = [
    (Rocky,       1.0,  15.0, 100.0,  0.0), // Flat beyond peak
    (Transitional, 2.0,  20.0, 150.0,  0.0),
    (Volatile,    5.0,  50.0, 500.0, -0.3), // Slight increase outward
    (Giant,       3.0, 1000.0, 5000.0, 0.3), // Bimodal: hot + cold
];
```

### Architecture Sampling by Spectral Type
```rust
fn sample_architecture(spectral: &str, metallicity: f64) -> SystemArchitecture {
    let giant_boost = 10.0_f64.powf(2.0 * metallicity);
    
    match spectral {
        "M" => weighted_choice(&[
            (CompactMulti, 0.45),
            (Sparse, 0.30),
            (Mixed, 0.20),
            (GiantDominated, 0.05 * giant_boost),
        ]),
        "K" => weighted_choice(&[
            (CompactMulti, 0.30),
            (Mixed, 0.30),
            (Sparse, 0.30),
            (GiantDominated, 0.10 * giant_boost),
        ]),
        "G" => weighted_choice(&[
            (CompactMulti, 0.25),
            (Mixed, 0.30),
            (Sparse, 0.30),
            (GiantDominated, 0.15 * giant_boost),
        ]),
        "F" | "A" => weighted_choice(&[
            (CompactMulti, 0.15),
            (Mixed, 0.25),
            (Sparse, 0.40),
            (GiantDominated, 0.20 * giant_boost),
        ]),
        _ => Sparse,
    }
}
```

---

## 9. Validation Targets

When generating 10,000+ systems, verify:

### Inner System (Kepler-constrained)
1. **Radius Distribution**: Bimodal with peaks at ~1.3 R⊕ and ~2.4 R⊕, gap at 1.5–2.0 R⊕
2. **Period Distribution**: Flat or rising to ~10d, then flat to ~100d, then gradual decline
3. **Giant-Metallicity**: Clear correlation, P(giant) should increase ~10× from [Fe/H]=-0.5 to +0.5
4. **M Dwarf Compactness**: M dwarf systems should average more planets, smaller radii
5. **Hot Jupiter Rate**: ~0.5–1% of FGK stars
6. **Hot Jupiter Loneliness**: Systems with hot Jupiters should have fewer inner companions

### Outer System (RV/DI/Microlensing-constrained)
7. **Cold Giant Rate**: ~10–20% of FGK stars should have a cold Jupiter (3–10 AU)
8. **Cold Giant Peak**: Giant occurrence should peak at 2–3× snow line, decline beyond
9. **Ice Giant Rate**: ~25–40% of systems should have at least one ice giant (10–40 AU)
10. **Ice Giant Multiplicity**: Of systems with ice giants: ~50% have 1, ~40% have 2, ~10% have 3
11. **Wide Companion Rate**: ~1–3% of systems should have a wide (50–300 AU) massive companion
12. **Metallicity for Ice Giants**: Weak or no correlation (unlike gas giants)
13. **Solar System Analogs**: ~5–10% of G-star systems should resemble Solar System architecture

---

## Key References

### Inner System / Kepler
1. **Fressin et al. (2013)** ApJ 766, 81 - False positive rates and occurrence
2. **Petigura et al. (2013)** PNAS 110, 19273 - Earth-size planet prevalence
3. **Howard et al. (2012)** ApJS 201, 15 - Planet occurrence within 0.25 AU
4. **Burke et al. (2015)** ApJ 809, 8 - Terrestrial planet occurrence rates
5. **Fulton et al. (2017)** AJ 154, 109 - Radius gap discovery
6. **Mulders et al. (2015, 2018)** - Stellar mass dependence, EPOS
7. **Hsu et al. (2019)** AJ 158, 109 - DR25 + Gaia occurrence rates
8. **Dressing & Charbonneau (2015)** ApJ 807, 45 - M dwarf occurrence
9. **Fischer & Valenti (2005)** ApJ 622, 1102 - Giant planet metallicity correlation
10. **Zhu & Dong (2021)** ARA&A 59, 291 - Architecture review

### Outer System / Radial Velocity
11. **Mayor et al. (2011)** arXiv:1109.2497 - HARPS survey, giant planet statistics
12. **Cumming et al. (2008)** PASP 120, 531 - Period-mass distribution from Keck
13. **Fernandes et al. (2019)** ApJ 874, 81 - Cold giant occurrence vs semi-major axis
14. **Wittenmyer et al. (2020)** MNRAS 492, 377 - Anglo-Australian Planet Search

### Microlensing
15. **Cassan et al. (2012)** Nature 481, 167 - One or more bound planets per star
16. **Suzuki et al. (2016)** ApJ 833, 145 - Mass ratio function beyond snow line
17. **Gould et al. (2010)** ApJ 720, 1073 - Microlensing planet statistics

### Direct Imaging
18. **Nielsen et al. (2019)** AJ 158, 13 - GPIES survey results
19. **Vigan et al. (2021)** A&A 651, A72 - SPHERE/SHINE demographics
20. **Bowler (2016)** PASP 128, 102001 - Direct imaging review
