# Load project environment
global path = dirname(@__FILE__)
import Pkg
Pkg.activate("$path")
Pkg.instantiate()

using CSV, DataFrames, Dates, DelimitedFiles, Distributed, Distributions, FileIO
using JLD2, MAT, Optim, OrderedCollections
using Printf, Random, Roots, SparseArrays, Statistics, UnPack
println("Packages loaded!")

## Build output folders if don't exist
global LOGS   = "$path/../output/log"
global OUTPUT = "$path/../output/data"
global INPUT  = "$path/../input"

!isdir("$path/../output/")	&& run(`mkdir $path/../output/`)
!isdir("$LOG/")				&& run(`mkdir $LOG/`)
!isdir("$OUTPUT/")			&& run(`mkdir $OUTPUT/`)
!isdir("$OUTPUT/../plots")	&& run(`mkdir $OUTPUT/../plots/`)
!isdir("$OUTPUT/../tables")	&& run(`mkdir $OUTPUT/../tables/`)
