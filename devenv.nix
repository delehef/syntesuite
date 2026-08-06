{ pkgs, lib, config, inputs, ... }:

{
  cachix.enable = false;
  
  packages = [ pkgs.git pkgs.sqlite pkgs.git-cliff ];

  languages.rust.enable = true;
}
