import type { Planet, PlanetClass } from "@magrathea/planetary-wasm";
import { match, P } from "ts-pattern";

/**
 * Get a color class for a planet based on its equilibrium temperature and class
 *
 * Returns Tailwind color class names for consistent, debuggable styling.
 *
 * Temperature ranges (Kelvin):
 * - Very Hot (>1000K): Orange/Red (lava worlds, hot Jupiters)
 * - Hot (400-1000K): Yellow/Orange (Venus-like, hot Neptunes)
 * - Temperate (200-400K): Blue/Green (Earth-like zone)
 * - Cold (50-200K): Cyan/Ice Blue (outer system, ice giants)
 * - Very Cold (<50K): Dark Blue (Kuiper belt objects)
 */
export function getPlanetColor(planet: Planet): string {
  const temp = planet.equilibriumTemp || 0;
  const planetClass = planet.class;

  // Class-specific color adjustments using Tailwind colors
  // Using lighter shades (200-400) for pale, pastel appearance
  return (
    match({ class: planetClass, temp })
      .with(
        { class: "Giant", temp: P.when((t) => t > 1000) },
        () => "bg-orange-200"
      )
      // Hot Jupiter
      .with(
        { class: "Giant", temp: P.when((t) => t > 400) },
        () => "bg-amber-200"
      )
      // Warm giant
      .with(
        { class: "Giant", temp: P.when((t) => t > 150) },
        () => "bg-amber-300"
      )
      // Cool giant (Jupiter-like)
      .with({ class: "Giant" }, () => "bg-orange-300") // Cold giant

      .with(
        { class: "Volatile", temp: P.when((t) => t > 400) },
        () => "bg-cyan-200"
      )
      // Hot Neptune
      .with(
        { class: "Volatile", temp: P.when((t) => t > 150) },
        () => "bg-sky-200"
      )
      // Warm Neptune
      .with(
        { class: "Volatile", temp: P.when((t) => t > 50) },
        () => "bg-blue-300"
      )
      // Ice giant (Neptune blue)
      .with({ class: "Volatile" }, () => "bg-blue-400")
      // Cold ice giant
      .with(
        { class: "Transitional", temp: P.when((t) => t > 600) },
        () => "bg-orange-200"
      )
      // Hot sub-Neptune
      .with(
        { class: "Transitional", temp: P.when((t) => t > 400) },
        () => "bg-orange-100"
      )
      // Warm
      .with(
        { class: "Transitional", temp: P.when((t) => t > 250) },
        () => "bg-sky-100"
      ) // Temperate (potential ocean)
      .with(
        { class: "Transitional", temp: P.when((t) => t > 150) },
        () => "bg-blue-200"
      )
      // Cool
      .with({ class: "Transitional" }, () => "bg-cyan-400")
      // Cold
      .with(
        { class: "Compact", temp: P.when((t) => t > 1000) },
        () => "bg-red-300"
      )
      // Lava world
      .with(
        { class: "Compact", temp: P.when((t) => t > 500) },
        () => "bg-red-200"
      ) // Hot rocky
      .with(
        { class: "Compact", temp: P.when((t) => t > 350) },
        () => "bg-pink-100"
      )
      // Warm rocky (Mars-like)
      .with(
        { class: "Compact", temp: P.when((t) => t > 280) },
        () => "bg-green-200"
      )
      // Temperate (Earth green)
      .with(
        { class: "Compact", temp: P.when((t) => t > 200) },
        () => "bg-slate-200"
      )
      // Cool
      .with(
        { class: "Compact", temp: P.when((t) => t > 100) },
        () => "bg-teal-200"
      )
      // Cold (icy surface)
      .with({ class: "Compact" }, () => "bg-slate-200")
      // Very cold (Pluto-like)
      .otherwise(() => "bg-gray-100")
  ); // Fallback gray
}

/**
 * Get a simple color based only on planet class (fallback if temperature unavailable)
 */
export function getPlanetColorByClass(planetClass: PlanetClass): string {
  return match(planetClass)
    .with("Giant", () => "bg-amber-300") // Tan (Jupiter-like)
    .with("Volatile", () => "bg-blue-300") // Neptune blue
    .with("Transitional", () => "bg-sky-200") // Sky blue
    .with("Compact", () => "bg-green-200") // Earth green
    .otherwise(() => "bg-gray-200");
}

/**
 * Get a temperature category label
 */
export function getTemperatureCategory(temp: number): string {
  if (temp > 1000) return "Very Hot";
  if (temp > 400) return "Hot";
  if (temp > 200) return "Temperate";
  if (temp > 50) return "Cold";
  return "Very Cold";
}

/**
 * Format temperature for display
 */
export function formatTemperature(temp: number): string {
  return `${Math.round(temp)} K`;
}
