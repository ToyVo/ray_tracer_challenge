import esbuild from "esbuild";
import { wasmLoader } from "esbuild-plugin-wasm";

await esbuild.build({
    entryPoints: ["browser/index.ts"],
    bundle: true,
    minify: true,
    sourcemap: false,
    format: "esm",
    outdir: "./public/dist",
    plugins: [wasmLoader()],
});
