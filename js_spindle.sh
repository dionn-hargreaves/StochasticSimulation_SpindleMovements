#!/bin/bash --login
#$ -cwd               # Application will run from current working directory
#$ -N julia_stochasticSpindle    # Name given to batch job (optional)
#$ -m bea
#$ -M dionn.hargreaves@postgrad.manchester.ac.uk



~/julia-1.7.2/bin/julia intro.jl




