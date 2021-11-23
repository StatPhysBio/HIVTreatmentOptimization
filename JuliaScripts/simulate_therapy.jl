include("file_access.jl")
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