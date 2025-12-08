import type { PlanetType } from "@magrathea/planetary-wasm";
import {
  scaleSemiMajorAxis,
  scaleUnifiedMass,
} from "../../../canvas/src/scale";
import { toEarthMasses } from "../utils/units";

type Props = {
  inclination: number;
  mass: number;
  maxWidth: number;
  planetType: PlanetType;
  semiMajorAxis: number;
};

export function PlanetCard({ mass, semiMajorAxis, maxWidth }: Props) {
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
  // planetType.type === "kuiperBeltObject" ? 8 * inclination - 16 : radius / 2;

  return (
    <div
      className="bg-white rounded-full absolute"
      style={{
        width: `${radius}px`,
        height: `${radius}px`,
        left: `${left}px`,
        top: `-${top}px`,
        background: "white",
      }}
    />
  );
}
