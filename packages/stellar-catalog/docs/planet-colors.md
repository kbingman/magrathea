# Planet Color Scheme

The planet visualization uses temperature-based coloring that varies by planet class for realistic appearance.

## Color Logic

Colors are derived from:
1. **Planet class** (Giant, Volatile, Transitional, Compact)
2. **Equilibrium temperature** in Kelvin
3. **Tailwind CSS color palette** for consistency
4. **80% opacity** for softer, more subtle appearance

This creates scientifically plausible colors - hot planets glow red/orange, temperate planets show blues/greens, and cold planets are dark blue/slate.

## Temperature Ranges

| Range | Category | Examples |
|-------|----------|----------|
| >1000K | Very Hot | Lava worlds, hot Jupiters |
| 400-1000K | Hot | Venus-like, hot Neptunes |
| 200-400K | Temperate | Earth-like zone |
| 50-200K | Cold | Outer system, ice giants |
| <50K | Very Cold | Kuiper belt objects |

## Colors by Class

All colors use Tailwind's palette with 80% opacity (`cc` suffix).

### Giant (Gas Giants)
- **>1000K**: `orange-500` - Orange (hot Jupiter)
- **400-1000K**: `amber-400` - Amber (warm giant)
- **150-400K**: `amber-600` - Dark amber (cool giant, Jupiter-like)
- **<150K**: `orange-800` - Dark orange (cold giant)

### Volatile (Ice Giants, Neptunes)
- **>400K**: `cyan-400` - Cyan (hot Neptune)
- **150-400K**: `sky-400` - Sky blue (warm Neptune)
- **50-150K**: `blue-500` - Blue (ice giant, Neptune)
- **<50K**: `blue-700` - Deep blue (cold ice giant)

### Transitional (Sub-Neptunes, Mini-Neptunes)
- **>600K**: `orange-400` - Orange (hot sub-Neptune)
- **400-600K**: `orange-300` - Light orange (warm)
- **250-400K**: `sky-300` - Light sky blue (potential ocean worlds)
- **150-250K**: `blue-400` - Blue (cool)
- **<150K**: `cyan-600` - Dark cyan (cold)

### Compact (Rocky/Terrestrial)
- **>1000K**: `red-600` - Red (lava world)
- **500-1000K**: `red-400` - Light red (hot rocky)
- **350-500K**: `pink-300` - Pink (Mars-like)
- **280-350K**: `green-500` - Green (habitable zone, Earth-like)
- **200-280K**: `slate-500` - Slate (cool)
- **100-200K**: `teal-500` - Teal (icy surface)
- **<100K**: `slate-700` - Dark slate (Pluto-like)

## Usage

```typescript
import { getPlanetColor } from '../utils/planet-colors';

const color = getPlanetColor(planet);
// Returns hex color like "#90be6d"
```

## Design Notes

- **Giant planets** use warmer earth tones (browns, tans) for cold giants, transitioning to oranges/yellows for hot Jupiters
- **Ice giants** use cyan/blue palette reflecting their methane atmospheres
- **Transitional planets** vary widely - could be mini-Neptunes (blues) or super-Earths with atmospheres
- **Rocky planets** show the most variation: lava worlds (red), desert worlds (pink/tan), Earth-like (green), icy (blue/teal)

The color scheme balances scientific plausibility with visual clarity.
