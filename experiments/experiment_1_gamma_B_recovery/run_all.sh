#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

Rscript code/01_fit_constant_B4.R
Rscript code/02_fit_piecewise_B_fixed_B0.R
Rscript code/03_fit_piecewise_B_estimate_B0_sigma.R
