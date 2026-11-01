#!env bash
set -o pipefail
set -e # Exit on error

# Rascal is not packaged on conda; install the 0.7.0 release from GitHub
tmp=$(mktemp -d)
curl -sfL -o "$tmp/rascal-0.7.0.tar.gz" https://github.com/crukci-bioinformatics/rascal/archive/refs/tags/0.7.0.tar.gz
R CMD INSTALL --no-docs "$tmp/rascal-0.7.0.tar.gz"
rm -rf "$tmp"
