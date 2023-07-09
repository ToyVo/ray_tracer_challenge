{
  description = "my project description";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.rust-overlay.url = "github:oxalica/rust-overlay";

  outputs = { self, nixpkgs, rust-overlay }:
    let
      inherit (nixpkgs) lib;
      forAllSystems = lib.genAttrs lib.systems.flakeExposed;
    in
    {
      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; overlays = [ rust-overlay.overlays.default ]; };
        in
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              cargo
              cargo-watch
              clippy
              rustfmt
              rust-analyzer
              wasm-pack
              nodejs
              trunk
              rust-bin.beta.latest.default
            ];
          };
        }
      );
    };
}
