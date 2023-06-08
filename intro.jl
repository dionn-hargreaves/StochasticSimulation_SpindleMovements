using DrWatson
@quickactivate "Stochastic-spindle-sim"



@info "Loading test parameters"
include("scripts/TestParameters.jl")

@info "Precompiling project"

push!(LOAD_PATH,"src")
using Revise
using BenchmarkTools
using StochasticSpindleSim

stochasticSpindleSim(Notes,NumGenerators,NumStates,finalTime,burnTime,maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)
