# *****************************
# Code for "Match Quality, Search, and the Internet Market for Used Books" (2021)
# by Glenn and Sara Fisher Ellison
# ...
# Julia implementation (Aug 2021 - June 2022) by Reca Sarfati (sarfati@mit.edu)
# Based on MATLAB code from Masao Fukui, Hongkai Zhang, [others]
# *****************************
using CSV, DataFrames, Distributions, FixedEffectModels, Plots, RegressionTables
using Roots, Statistics

### Build output folders if don't exist
!isdir("plots")  && run(`mkdir plots/`)
!isdir("tables") && run(`mkdir tables/`)

### Users should Ctrl-F for all TODO's in main.jl for personal use.

### Acquire parallel workers (specific to MIT server)
### TODO: Accomodate code for personal computing environment.



### 
