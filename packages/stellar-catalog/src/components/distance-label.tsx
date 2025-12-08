import type { PropsWithChildren } from "react";
import { scaleSemiMajorAxis } from "@magrathea/canvas";

const OFFSET_LEFT = 60;

type Props = {
  semiMajorAxis: number;
  maxOutput: number;
};

export function DistanceLabel({
  children,
  semiMajorAxis,
  maxOutput,
}: PropsWithChildren<Props>) {
  const left = Math.round(
    scaleSemiMajorAxis(semiMajorAxis, { maxOutput }) + OFFSET_LEFT
  );

  return (
    <div
      className="text-xs text-neutral-500 absolute top-0"
      style={{
        left: `${left}px`,
      }}
    >
      {children}
    </div>
  );
}
