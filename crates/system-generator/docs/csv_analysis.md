# Planet Generation: CSV Analysis (1000 Systems, 7073 Planets)

## Mass-Radius Relationships

The scatter fix appears to have worked. Densities are now physically plausible:

| Category | Mean Density (Earth=1) | Expected | Status |
|----------|------------------------|----------|--------|
| Rocky planets | 0.85 | 0.7–1.5 | ✓ Good |
| Transitional | 0.58 | 0.3–0.7 | ✓ Good |
| Ice Giants | 0.63 | ~0.3 (Neptune) | Slightly high |
| Gas Giants | 0.21 | ~0.24 (Jupiter) | ✓ Good |

Sample Terran planets show sensible values (M=1.13, R=1.10 → ρ=0.84).

## Giant Planet Occurrence by Stellar Mass

| Stellar Type | n | Gas Giant | Ice Giant | Expected Gas | Status |
|--------------|---|-----------|-----------|--------------|--------|
| Late M (M6–M9) | 295 | 0.7% | 13.9% | <1–2% | ✓ |
| Mid M (M3–M5) | 313 | 4.5% | 18.8% | 3–5% | ✓ |
| Early M (M0–M2) | 170 | 14.1% | 27.1% | 3–5% | ⚠️ High |
| K dwarfs | 99 | 16.2% | 34.3% | ~10% | ⚠️ Bit high |
| G dwarfs | 44 | 29.5% | 31.8% | 10–15% | ⚠️ High |
| F/A dwarfs | 34 | 38.2% | 41.2% | 15–20% | ⚠️ High |

Late M suppression is working well. Early M and higher-mass stars are forming giants too readily.

## Concerning Cases

14 gas giants formed around stars below 0.3 M☉, including:

- M8 (0.088 M☉) → Gas Giant 295 M⊕ — planet-to-star mass ratio of 1:100
- M5 (0.154 M☉) → Gas Giant 402 M⊕ — ratio of 1:128

These are physically implausible. A gas giant at 1:100 mass ratio around an M8 would likely have disrupted disk dynamics before forming.

## Ice Giant Rates

Ice giants appear too frequently around M-dwarfs:

| Stellar Type | Ice Giant Rate | Expected |
|--------------|----------------|----------|
| Late M | 13.9% | ~5% |
| Mid M | 18.8% | ~10% |
| Early M | 27.1% | ~15% |

The suppression may be applied to gas giants but not ice giants, or the ice giant threshold is too weak.

## Planet Type Distribution

The type diversity looks good:

| Type | Count | Mass Range (M⊕) | Radius Range (R⊕) |
|------|-------|-----------------|-------------------|
| KBO | 1970 | <0.01 | 0.05–0.10 |
| Dwarf Planet | 1436 | 0.0001–0.01 | 0.09–0.35 |
| Sub-Earth | 1360 | 0.01–0.50 | 0.26–0.92 |
| Mini-Neptune | 455 | 2.0–5.0 | 1.2–2.4 |
| Desert | 445 | 0.5–2.0 | 0.8–1.4 |
| Frozen | 331 | <2.0 | 0.05–1.4 |
| Terran | 264 | 0.5–2.0 | 0.8–1.3 |
| Ice Giant | 233 | 5–158 | 2.1–6.0 |
| Gas Giant | 81 | 161–2366 | 10–22 |

Rare types (Hot Jupiter, Puffy Saturn, Oceanic) appear in small numbers as expected.

## Recommendations

1. **Strengthen early M suppression** — gas giants at 14% for M0–M2 is too high
2. **Add ice giant suppression** — currently only gas giants seem suppressed
3. **Cap planet-to-star mass ratio** — prevent 1:100 cases around very low mass stars
4. **Review G/K rates** — 30% gas giant occurrence is roughly 2–3x too high
