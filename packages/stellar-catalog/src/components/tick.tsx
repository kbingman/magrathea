import { scaleSemiMajorAxis } from "@magrathea/canvas";

export function HalfTicks({
  scale,
  maxOutput,
}: {
  scale: number;
  maxOutput: number;
}) {
  const step = scale / 10;
  return (
    <>
      {Array.from({ length: 8 }).map((_, i) => (
        <div
          key={`tick-${step * (i + 2)}`}
          className="absolute top-0"
          style={{
            left: `${scaleSemiMajorAxis(step * (i + 2), {
              maxOutput,
            })}px`,
          }}
        >
          <div className="absolute h-1 -top-0.5 w-px bg-neutral-700 text-xs" />
        </div>
      ))}
    </>
  );
}

type Props = {
  semiMajorAxis: number;
  maxOutput: number;
};

export function Tick({ semiMajorAxis, maxOutput }: Props) {
  const left = scaleSemiMajorAxis(semiMajorAxis, {
    maxOutput,
  });

  return (
    <>
      <HalfTicks scale={semiMajorAxis} maxOutput={maxOutput} />
      <div
        className="absolute top-0"
        style={{
          left: `${left}px`,
        }}
      >
        <div className="absolute h-2 -top-1 w-px bg-neutral-700" />
      </div>
    </>
  );
}
