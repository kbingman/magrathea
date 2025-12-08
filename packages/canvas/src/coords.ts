/**
 * @param ctx Canvas Context
 * @param scale the current scale
 * @param coord the X coordinate
 * @returns the translated coordinate
 */
export function translateX(
  ctx: CanvasRenderingContext2D,
  scale: number,
  coord: number
) {
  return coord * scale + ctx.canvas.width / 2;
}

/**
 * @param ctx Canvas Context
 * @param scale the current scale
 * @param coord the X coordinate
 * @returns the translated coordinate
 */
export function translateY(
  ctx: CanvasRenderingContext2D,
  scale: number,
  coord: number
) {
  return ctx.canvas.height / 2 - coord * scale;
}
