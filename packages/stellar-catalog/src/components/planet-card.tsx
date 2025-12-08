import type { Planet } from "@magrathea/planetary-wasm";
import {
  scaleSemiMajorAxis,
  scaleUnifiedMass,
} from "../../../canvas/src/scale";
import { toEarthMasses } from "../utils/units";
import { getPlanetColor } from "../utils/planet-colors";

type Props = {
  planet: Planet;
  maxWidth: number;
};

export function PlanetCard({ planet, maxWidth }: Props) {
  const { mass, semiMajorAxis, class: planetClass, planetType } = planet;

  const left = Math.round(
    scaleSemiMajorAxis(semiMajorAxis, {
      maxOutput: maxWidth,
    })
  );
  const radius = Math.round(
    scaleUnifiedMass(toEarthMasses(mass), {
      maxOutput: 8,
      minInput: 0.1,
      minOutput: 2,
    })
  );

  const top = radius / 2;
  const colorClass = getPlanetColor(planet);

  return (
    <div
      className={`rounded-full absolute ${colorClass}`}
      title={`${planetClass} - ${planetType.name}`}
      style={{
        width: `${radius}px`,
        height: `${radius}px`,
        left: `${left}px`,
        top: `-${top}px`,
      }}
    />
  );
}
