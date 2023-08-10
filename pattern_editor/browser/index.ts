if (DEBUG) {
    new EventSource('/esbuild').addEventListener('change', () => location.reload())
}
import { Universe } from "../pkg/pattern_editor";
import { memory } from "../pkg/pattern_editor_bg.wasm";
const canvas = document.getElementById("canvas") as HTMLCanvasElement
const universe = Universe.new();
canvas.height = universe.canvas_height();
canvas.width = universe.canvas_width();
const ctx = canvas.getContext("2d") as CanvasRenderingContext2D;
const renderLoop = () => {
    const dataPtr = universe.data();
    const data = new Uint8ClampedArray(memory.buffer, dataPtr, universe.size() * 4);
    const imageData = new ImageData(data, universe.canvas_width(), universe.canvas_height());
    ctx.putImageData(imageData, 0, 0);
    universe.tick();
    requestAnimationFrame(renderLoop);
};
requestAnimationFrame(renderLoop);
