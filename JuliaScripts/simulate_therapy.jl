include("file_access.jl")
# This is a list of scripts related to the clinical trial data 
# and creating histograms realted to three trials
# We also simulate trial related quantities, such as the role of mutations
# and the traces of the  viral fractions.
##
simulate_trial(trial; kwds...) = simulate_trial_bayes(trial,trial_antibodies(trial);kwds...)
histogramfile = datapath*"/histogramdata.h5"

## Write the transformed histograms out
h5open(histogramfile,"w") do fid
    for trial in trial_list
        fid[trial] = simulate_trial(trial;diversity_multiplier = 2.07)
    end
end
# write the untransformed hisotgrams out
h5open(histogramfile,"r+") do fid
    for trial in trial_list
        fid[trial*"ref"] = simulate_trial(trial;diversity_multiplier = 1.0)
    end
end


# Effect of mutations

function mut_compare_bayes(trial) 
    profile = ab_profile_list_bayes(trial_antibodies(trial)) # share parameter disorder
    m_true = trial_rebound_prob(get_start_theta(trial) .* reservoir_factor, profile;
        n_samples = 100, mutations = true)
    m_false = trial_rebound_prob(get_start_theta(trial) .* reservoir_factor, profile;
        n_samples = 100, mutations = false)
    [m_true,m_false]
end

mut_compare(trial) = mean(mut_compare_bayes(trial) for _ in 1:20)

rebound_prob_1074 = mut_compare("10-1074")

rebound_prob_3BNC = mut_compare("3BNC")

rebound_prob_combo = mut_compare("combo")
##
fid = h5open(datapath*"/trialsimulations_mutvsnomut.h5", "w")
fid["10-1074"] =  rebound_prob_1074
fid["3BNC117"] = rebound_prob_3BNC
fid["combo"] =  rebound_prob_combo
close(fid)

## Simulaitons for visualization


function trial_traces(trial, ab_list; samples_per_patient = 10, 
	# save full traces for plotting trajectory swarm
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = 2,
	timepoints = 0:.2:56
	)
	theta_vec = get_start_theta(trial)
	trace = []
	out_mat = permutedims(hcat(collect(timepoints),collect(timepoints)))
	for (ii,θ_0) in enumerate(theta_vec .* diversity_multiplier)
		ab_profile_list = ab_profile_list_bayes(ab_list; selection_multiplier)
		vp = initialize_viral_population(θ_0, ab_profile_list; λ = 20.0)
        # λ = 20 increase the fineness of population discretization.
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end
		for jj in  1:samples_per_patient
			restart_viral_population!(vp)
			tr = virus_time_trace(vp, timepoints; breakpoint)
			trace =  hcat(map(tr) do x 
				[x...] ./ vp.capacity
			end...)
			out_mat = cat(out_mat, trace, dims = 3)
		end
	end
	return out_mat
end

tt = trial_traces("10-1074",["10-1074"], samples_per_patient = 100, diversity_multiplier=2.07)
fid = h5open(datapath*"/trialsimulations_traces.h5", "w")
fid["10-1074"] = tt
close(fid)