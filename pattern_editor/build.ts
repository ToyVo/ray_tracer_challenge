import esbuild from "esbuild";
import wasmpack from "esbuild-plugin-wasm-pack";

esbuild.build({
    entryPoints: ["browser/app.tsx"],
    bundle: true,
    minify: true,
    sourcemap: true,
    target: ["chrome58", "firefox57", "safari11", "edge16"],
    outdir: "./public/dist",
    plugins: [
        wasmpack({
            path: ".",
        }),
    ],
});
