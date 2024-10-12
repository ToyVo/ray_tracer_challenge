{
  description = "discord_bot, A Rust web server including a NixOS module";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-parts = {
      url = "github:hercules-ci/flake-parts";
      inputs.nixpkgs-lib.follows = "nixpkgs";
    };
    crate2nix.url = "github:nix-community/crate2nix";

    # Development

    devshell = {
      url = "github:numtide/devshell";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  nixConfig = {
    extra-trusted-public-keys = "eigenvalue.cachix.org-1:ykerQDDa55PGxU25CETy9wF6uVDpadGGXYrFNJA3TUs=";
    extra-substituters = "https://eigenvalue.cachix.org";
    allow-import-from-derivation = true;
  };

  outputs =
    inputs@{
      self,
      nixpkgs,
      flake-parts,
      crate2nix,
      devshell,
      ...
    }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      imports = [
        # devshell.flakeModule
        flake-parts.flakeModules.easyOverlay
      ];

      flake = {
      };

      perSystem =
        {
          self',
          system,
          pkgs,
          lib,
          config,
          ...
        }:
        let
          cargoToml = builtins.fromTOML (builtins.readFile ./ray_tracer/Cargo.toml);
          rev = self.shortRev or self.dirtyShortRev or "dirty";
          generatedCargoNix = crate2nix.tools.${system}.generatedCargoNix {
            name = "ray_tracer";
            src = ./.;
          };
          cargoNix = pkgs.callPackage "${generatedCargoNix}/default.nix" {
            buildRustCrateForPkgs =
              pkgs:
              pkgs.buildRustCrate.override {
                defaultCrateOverrides = pkgs.defaultCrateOverrides // {
                  ray_tracer = attrs: {
                    version = "${cargoToml.package.version}-${rev}";
                    buildInputs =
                      with pkgs.darwin.apple_sdk.frameworks;
                      lib.optionals pkgs.stdenv.isDarwin [
                        SystemConfiguration
                        CoreServices
                      ];
                    nativeBuildInputs = with pkgs; [
                      libiconv
                      pkg-config
                    ];
                  };
                };
              };
          };
        in
        {
          packages = {
            ray_tracer = cargoNix.workspaceMembers.ray_tracer.build;
            default = self'.packages.ray_tracer;
          };
          overlayAttrs = {
            inherit (self'.packages) ray_tracer;
          };
          devShells.default =
            pkgs.mkShell {
              shellHook = ''
                export RUST_LOG="ray_tracer=trace"
                export RUST_SRC_PATH=${pkgs.rustPlatform.rustLibSrc}
              '';
              buildInputs = with pkgs.darwin.apple_sdk.frameworks; [ SystemConfiguration ];
              nativeBuildInputs = with pkgs; [
                rustc
                pkg-config
                rustPlatform.bindgenHook
                libiconv
                cargo-watch
                systemfd
              ];
            };
        };
    };
}
