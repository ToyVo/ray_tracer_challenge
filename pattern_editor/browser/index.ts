if (DEBUG) {
    new EventSource("/esbuild").addEventListener("change", () =>
        location.reload(),
    );
}
import { PatternImage } from "../pkg/pattern_editor";
import { memory } from "../pkg/pattern_editor_bg.wasm";
const canvas = document.getElementById("canvas") as HTMLCanvasElement;
const universe = PatternImage.new();
const width = universe.width();
const height = universe.height();
const bufferLength = width * height * 4;
canvas.height = height;
canvas.width = width;
const ctx = canvas.getContext("2d") as CanvasRenderingContext2D;
const dataPtr = universe.data();
const data = new Uint8ClampedArray(memory.buffer, dataPtr, bufferLength);
const imageData = new ImageData(data, width, height);
ctx.putImageData(imageData, 0, 0);
