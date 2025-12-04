# Statistical Planet Generator Assessment

Analysis of 1000 systems each: IMF sample (Kroupa distribution, M-dwarf dominated) and Solar Analog sample (G-type only).

## Overview Comparison

|                              | IMF Sample       | Solar Analog |
|------------------------------|------------------|--------------|
| Total planets                | 3,607            | 2,817        |
| Planets/system               | 4.13             | 3.40         |
| Unique fingerprints          | 16.5%            | 16.8%        |

### Planet Type Distribution

| Type          | IMF    | Solar  | Notes |
|---------------|--------|--------|-------|
| Sub-Earth     | 38.0%  | 32.2%  | High but realistic |
| Gas Giant     | 5.1%   | 8.3%   | |
| Ice Giant     | 12.3%  | 16.7%  | |
| Hot Jupiter   | 0.2%   | 0.4%   | Good - should be rare |

### Architecture

| Feature                    | IMF    | Solar  |
|----------------------------|--------|--------|
| Compact (≥4p, <1 AU)       | 34.4%  | 0.7%   |
| Extended (>5 AU)           | 18.8%  | 45.5%  |
| Multi-giant                | 2.3%   | 5.6%   |
| High eccentricity (e>0.4)  | 3.1%   | 3.0%   |

## Issues Identified

### 1. Giant Formation Around M-Dwarfs [CRITICAL]

**Problem:** 15-16% of M-dwarfs have gas giants, with planet/star mass ratios up to 15%. A 12 M_J planet around a 0.09 M☉ star is effectively a brown dwarf companion, not a planet.

**Literature:** Johnson et al. (2010) found giant occurrence scales as ~M_star^1.0. For M-dwarfs, giant occurrence should be ~2-5%, not 15%.

**Fix:**
```rust
// Scale giant probability with stellar mass
let giant_prob_base = /* from metallicity */;
let stellar_scaling = (stellar_mass / 1.0).powf(1.5); // α ~ 1-2
let giant_prob = giant_prob_base * stellar_scaling;

// Also cap planet mass relative to stellar mass
let max_planet_mass = stellar_mass * 317.8 * 0.01; // 1% of stellar mass
```

### 2. Volatile Class Too Broad

**Problem:** "Ice Giants" span 5.8 - 160 M⊕. The upper end (160 M⊕ = 0.5 M_J) is Saturn-class, not Neptune-class.

**Fix:** Split the Volatile class:
- Ice Giant: 10-30 M⊕ (Neptune-class)
- Sub-Saturn: 30-95 M⊕ 
- Keep Saturn-class (95-160 M⊕) distinct from Jupiter-class

### 3. Sub-Earth Prevalence (32-38%)

**Observation:** Many systems look like "5 Sub-Earths + 1 mini-Neptune". This is probably realistic (Kepler can't detect most sub-Earths) but contributes to fingerprint clustering.

**Consider:** Minimum mass threshold? Or just accept this as "debris" population that adds realism but not interest.

### 4. Fingerprint Clustering

**Observation:** Only 16-17% unique fingerprints. Top 5 account for ~25% of systems.

**Mitigation:** Within-fingerprint diversity is real - orbital spacing varies 1.5× to 28× within same fingerprint. The problem is less severe than the number suggests.

## What's Working Well

- ✓ **Metallicity-giant correlation:** 1% at [Fe/H]<-0.3 → 49% at [Fe/H]>+0.3
- ✓ **Hot Jupiter placement:** All <0.1 AU, T>1500K, ~0.2-0.4% occurrence
- ✓ **Temperature-type consistency:** Lava Worlds hot, Ice Giants cold
- ✓ **Wide separation companions:** 1-2%, all gas giants (matches direct imaging)
- ✓ **M-dwarf orbital structure:** Median 0.06 AU, concentrated near star
- ✓ **Solar analog outer systems:** 46% have planets >5 AU
- ✓ **Resonant pairs:** ~18% of systems have near-resonant pairs
- ✓ **Eccentricity outliers:** 3% with e>0.4

## Notable Systems Found

### Solar System Analogues (3.3%)
Systems with HZ terrestrial + Jupiter-zone giant + distant ice giant.

### High Eccentricity (3%)
Giant planets with e > 0.4, suggesting past dynamical interactions.

### Hot Jupiters with Companions (rare)
Only 2 systems - consistent with observation that hot Jupiters tend to be lonely.

### Multi-Giant Systems
5.6% of solar analogs have 2+ giants - provides architectural diversity.

## Recommended Priority Fixes

1. **[HIGH]** Scale giant formation with stellar mass
2. **[MEDIUM]** Cap planet mass at ~1% stellar mass  
3. **[LOW]** Split Volatile class for better mass resolution
4. **[CONSIDER]** Sub-Earth threshold for "interesting" systems

## The "Boring" Verdict

**Not as bad as feared.** While fingerprint clustering exists, the orbital architectures within fingerprints do vary. Systems that "look the same" on paper actually differ in spacing, mass distribution, and resonance structure.

The real diversity killers are:
- Sub-Earth dominance (every system has several)
- Single-planet Sparse systems (127 in IMF sample)
- Lack of dramatic outliers (migration scars, instability survivors)

For storytelling/worldbuilding purposes, these systems are serviceable. For scientific comparison, they need the giant-stellar-mass fix.
