using HDF5: eachindex




include("file_access.jl")

## 
div_mult = 2 .^(-1:.05:2.5)
sel_mult = 10 .^(-1:.05:1)

## Current best data for reservoir
mult_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = 1.0, diversity_multiplier = d, n_pars = 100)
			for d in div_mult) 
			for trial in ["10-1074","3BNC","combo"])
##
bestind = argmin(ind -> sum(mult_out, dims = 2)[ind], eachindex(div_mult))
minind = argmin(ii -> log(div_mult[ii])^2, eachindex(div_mult))
stat_diversity = sum(mult_out, dims = 2)[minind] - sum(mult_out, dims = 2)[bestind] # difference between min and optimal
best_diversity = div_mult[bestind] # actually optimal diversity
##
sel_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = s, diversity_multiplier = best_diversity, n_pars = 100)
			for s in sel_mult) 
			for trial in ["10-1074","3BNC","combo"])
##


##

h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_selection.h5", "w") do fid
	fid["10-1074"] = out[:,1]
	fid["3BNC117"] = out[:,2]
	fid["combo"] = out[:,3]
	fid["sel_mult"] = sel_mult
end
