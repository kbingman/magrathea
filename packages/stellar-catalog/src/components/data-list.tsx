import type { ReactNode } from "react";

type Props = {
  data: [ReactNode, ReactNode][];
  title: string;
};

export function DataList({ data, title }: Props) {
  return (
    <div className="p-2 border border-amber-900 my-2 grid gap-1.5">
      <h2 className="text-sm uppercase">{title}</h2>
      <dl className="text-xs uppercase grid grid-cols-2 gap-1 justify-between">
        {data.map(([k, v]) => (
          <>
            <dt className="text-left">{k}</dt>
            <dd className="text-right">{v}</dd>
          </>
        ))}
      </dl>
    </div>
  );
}
