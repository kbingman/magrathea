import type {
  MainSequenceStar,
  Planet,
  StellarObject,
} from "@magrathea/planetary-wasm";
import { useState } from "react";
import { scaleUnifiedMass } from "../../../canvas/src/scale";
import { toEarthMasses } from "../utils/units";
import { PlanetCard } from "./planet-card";
import { Tick } from "./tick";

type Props = {
  id: string;
  name: string;
  stars: StellarObject[];
  planets: Planet[];
  maxWidth: number;
};

export function SystemCard({ id, name, planets, stars, maxWidth }: Props) {
  const [showInfo, _setShowInfo] = useState(false);
  const { color, luminosityClass, mass, spectralType, subtype } =
    stars[0] as MainSequenceStar;

  const stellarRadius = Math.round(
    scaleUnifiedMass(mass, {
      maxOutput: 32,
      minInput: 0.05,
      minOutput: 4,
    }),
  );
  const stellarColor = `rgb(${color.r}, ${color.g}, ${color.b})`;

  return (
    <div id={`system-${id}`} className="flex items-center gap-6 relative">
      <div className="text-xs whitespace-nowrap">{name}</div>
      <div className="relative my-2 h-px bg-neutral-800 w-full">
        <Tick semiMajorAxis={1} maxOutput={maxWidth} />
        <Tick semiMajorAxis={10} maxOutput={maxWidth} />
        <Tick semiMajorAxis={100} maxOutput={maxWidth} />
        <Tick semiMajorAxis={1000} maxOutput={maxWidth} />
        <div
          className="bg-white rounded-full absolute"
          style={{
            width: `${stellarRadius}px`,
            height: `${stellarRadius}px`,
            left: `-${stellarRadius / 2}px`,
            top: `-${stellarRadius / 2}px`,
            background: stellarColor,
          }}
        />
        {planets.map((planet) => (
          <PlanetCard
            key={`icon-${planet.id}`}
            planet={planet}
            maxWidth={maxWidth}
          />
        ))}
      </div>
      {showInfo && (
        <div className="absolute -top-2 left-24 bg-neutral-950 text-xs p-2 border-amber-900 border z-10 w-48">
          <h2>
            {spectralType}
            {subtype}
            {luminosityClass} {mass.toFixed(3)}M☉
          </h2>
          <div className="pl-2">
            {planets
              .filter(
                (p) =>
                  p.planetType.type !== "kuiperBeltObject" &&
                  p.planetType.type !== "dwarfPlanet",
              )
              .map(({ id, planetType, mass }) => {
                console.log({ id });
                return (
                  <div key={id}>
                    {planetType.name} {toEarthMasses(mass).toFixed(3)}M⊕
                  </div>
                );
              })}
          </div>
        </div>
      )}
    </div>
  );
}
