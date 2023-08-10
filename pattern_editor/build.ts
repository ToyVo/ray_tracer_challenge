import esbuild from "esbuild";
import { wasmLoader } from "esbuild-plugin-wasm";
import wasmpack from "esbuild-plugin-wasm-pack";

esbuild.build({
    entryPoints: ["browser/index.ts"],
    bundle: true,
    minify: true,
    sourcemap: false,
    format: "esm",
    outdir: "./public/dist",
    plugins: [
        wasmLoader(),
        wasmpack({
            path: ".",
        }),
    ],
    define: {
        DEBUG: "false",
    },
});
