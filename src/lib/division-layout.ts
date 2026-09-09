// The text occupies 66% of the cell width. Keep its height (including entrance
// movement) within 64% of the cell height so all four corners fit inside the ellipse.
export function divisionCellHeight(width: number, textHeight: number): number {
  return Math.ceil(Math.max(width, (textHeight + 52) / 0.64));
}
