import test from 'node:test';
import assert from 'node:assert/strict';
import { divisionCellHeight } from '../src/lib/division-layout.ts';

test('narrow-screen copy and its entrance movement stay inside the final cell ellipse', () => {
  for (const width of [256, 288, 343, 398, 480]) {
    for (const textHeight of [160, 240, 320, 480, 640]) {
      const cellHeight = divisionCellHeight(width, textHeight);
      // Radius is 158 in a 320-unit SVG cell; test the furthest text corner
      // after its maximum 26px entrance translation.
      const radiusX = width * 158 / 320;
      const radiusY = cellHeight * 158 / 320;
      const cornerX = width * 0.66 / 2;
      const cornerY = textHeight / 2 + 26;
      assert.ok((cornerX / radiusX) ** 2 + (cornerY / radiusY) ** 2 < 1,
        `copy ${width} × ${textHeight} must remain within its cell`);
    }
  }
});
