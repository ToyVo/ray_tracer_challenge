import esbuild from "esbuild";
import { wasmLoader } from "esbuild-plugin-wasm";
import wasmpack from "esbuild-plugin-wasm-pack";

(async () => {
    const ctx = await esbuild.context({
        entryPoints: ["browser/index.ts"],
        bundle: true,
        minify: true,
        sourcemap: true,
        format: "esm",
        outdir: "./public/dist",
        plugins: [
            wasmLoader(),
            wasmpack({
                path: ".",
            }),
        ],
        define: {
            DEBUG: "true",
        },
    });

    await ctx.watch();

    const { host, port } = await ctx.serve({
        servedir: "./public",
    });

    console.log(`Server listening on http://${host}:${port}`);
})();
