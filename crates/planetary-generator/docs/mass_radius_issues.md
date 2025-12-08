# Mass-Radius Relationship Issues

## Core Problem

Mass and radius appear to be sampled independently rather than derived from a physically-motivated relationship. This produces planets with impossible densities for their assigned composition class.

## Specific Issues

### 1. Rocky Planets with Gas Giant Densities

213 rocky planets have R > 2.0 R⊕, with the worst cases showing:

| Type | Mass (M⊕) | Radius (R⊕) | Density (ρ⊕) | Problem |
|------|-----------|-------------|--------------|---------|
| Lava World | 1.85 | 3.72 | 0.036 | Less dense than Saturn |
| Desert | 1.81 | 3.65 | 0.037 | Requires massive H/He envelope |
| Terran | 1.77 | 3.30 | 0.049 | Cannot be rocky at this density |

A 2 M⊕ rocky planet should have R ≈ 1.4-1.5 R⊕, not 3+ R⊕.

### 2. Iron Worlds Not Actually Dense

Iron worlds should be the densest planetary class (ρ ≈ 1.5-2.5 ρ⊕), but the data shows:

- Density range: 0.04 to 22.0 ρ⊕
- Many iron worlds are fluffier than Earth
- Type assignment appears disconnected from actual mass-radius values

### 3. No Consistent M-R Relation Within Types

Terran planets show radius varying from 0.33× to 2.16× the expected value for their mass. This scatter is far too large to represent natural variation — it indicates radius is not conditioned on mass.

### 4. High-Density Small Planets

870 planets with M < 2 M⊕ have density > 2.0 ρ⊕. While some high-density planets exist (Mercury-like iron-rich bodies), this is too many. Most are typed as Sub-Earth, Desert, or Terran rather than Iron World.

## Expected Mass-Radius Relations

For reference, empirical and theoretical M-R relations:

**Rocky (Earth-like composition):**
```
M < 1 M⊕:   R = M^0.27
M > 1 M⊕:   R = M^0.55  (compression becomes significant)
```

**Composition variants:**
```
Iron-rich:    R = 0.80 × R_rocky  (Mercury-like)
Water-rich:   R = 1.15 × R_rocky  (lower bulk density)  
```

**Volatile-dominated (Neptunes, Jupiters):**
```
Different regime — radius depends on envelope mass fraction,
equilibrium temperature, and age (contraction)
```

## Recommended Fix

Radius should be derived from mass, not sampled independently:

1. Assign mass based on formation/type logic (current approach may be fine)
2. Compute baseline radius from mass using composition-appropriate M-R relation
3. Add modest scatter (σ ≈ 10-15%) to represent natural variation
4. For volatile-rich planets, account for envelope inflation from stellar heating

This ensures density is always physically plausible for the assigned composition class.

## Validation Check

After fixing, verify:
- All rocky planets have ρ > 0.7 ρ⊕
- Iron worlds have ρ > 1.3 ρ⊕
- No rocky planets exceed ~1.8 R⊕ (the radius gap)
- Transitional planets span the 1.8-3.5 R⊕ range appropriately
