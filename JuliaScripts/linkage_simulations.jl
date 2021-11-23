#using Pkg
#Pkg.add(url ="https://github.com/StatPhysBio/Tomoko")


using Distributed # older code that uses distributed rather than Tomoko.
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Tomoko
using Tomoko
using CSV
using Printf
using DataFrames


path = datapath = joinpath(@__DIR__,"../H5Output")
cd(path)

function save_pop_fit(path, result, par; name_fields = [:κ,:ν,:χ], p = 0)
    s = ["$(name)=$(@sprintf("%.0E", getfield(par,name)))_" for name in name_fields]
    CSV.write(string(path,s...,"pflip=$(@sprintf("%.0E", p))",".csv"),result)
end

## Observables 
# θ = 2*μ*Ne  = 2*(ν *λ / 2) κ/λ = ν κ
# rs = Δ/μ  = Δ/(ν *λ / 2)  =>  Δ = (ν *λ / 2) rs

global loci = 2^9
global lambda = 4
global kappa = 1000
for x in [0,0.1], p in [0,0.01], θ_ in 0.02
    β1 = zeros(loci)
    nu = θ_ / kappa # the first number sets the diversity 
    fixed_sites = 2^3:2^3:2^9 # 64 fixed sites running from Δ/μ = 8 to 512 in increments of 8
    for ii in fixed_sites
        β1[ii] += ii * (nu * lambda) /2
    end
    variable_sites = 2^4-1:2^4:2^9 # 32 variable sites
    for ii in variable_sites
        β1[ii] += .02 # 2% fitness difference (max - mean < 4 to avoid blow up)
    end
    true_fit =β1[fixed_sites]/((nu * lambda) /2)
    result = DataFrame(f = true_fit)
    par = PopRates(
        κ= kappa,
        χ= x, 
        ρ= 0.1, 
        ν= nu , # first number determines D
        loci = loci, 
        β1 = β1)
    for ii in 1:80
        dglist = pmap(x->run_sim(par, 0:10:2000 ; var_sites = variable_sites, flip_prob = p/length(variable_sites)), 1:11)
        fitness = pmap(ii->estimate_rs(dglist, locus = ii, timepoints = 1000:100:2000) , fixed_sites)
        result[Symbol("run",ii)]=fitness
    end
    save_pop_fit(path, result, par, p = p)
end
