
using Combinatorics
include("file_access.jl")

# trialsimulations_statistics.h5 is
# If available:
#  h5open(datapath*"trialsimulations_statistics.h5", "r") do fid
#      global reservoir_factor = read(fid["best_fit_multplier"])
#  end
##
global reservoir_factor = 2.07
## Combinations



function rebound_dict(ab_list, n_antibodies; quartile_estimator_samples = 4) 
    #generate a dictionary of therapies to exhaustive combinations 
    # of bayesian realizations of rebound probabilities
    out = Dict{Array{String,1},Array{Float64,1}}() 
    for combo in combinations( ab_list, n_antibodies)
        list = Float64[]
        for _ in 1:quartile_estimator_samples
			p = trial_rebound_prob(get_start_theta("all") .* reservoir_factor, 
                ab_profile_list_bayes(combo); n_samples = 50)
			push!(list, p)
        end
        out[combo] = list
        print(combo)
    end
    return out
end

function make_refinements!(dictionary; 
    top = 10, 
    addenda= [["10-1074"],["3BNC117"],["10-1074","3BNC117"]], # anything we want to refine for sure
    quartile_estimator_samples = 4)
    # Identify the top ten treatment candidates to construct more
    # simulations for
    relevant_addenda = intersect(collect(keys(dictionary)),addenda)
    k_list = sort(collect(keys(dictionary)), by = x->median(dictionary[x]))[1:min(end,top)]
    k_list = union(k_list,relevant_addenda) # make sure special cases are added
    for k in k_list
        for _ in 1:quartile_estimator_samples
            p = trial_rebound_prob(get_start_theta("all") .* reservoir_factor, 
                ab_profile_list_bayes(k); n_samples = 50)
            push!(dictionary[k], p)
        end
    end
end 

function write_rebound_h5(outdict)
    ii = length(first(keys(outdict)))
    h5open(datapath*"/therapy_ranking.h5","r+") do fid
        loc = create_group(fid, "$(ii)_antibodies")
        for k in keys(outdict)
            kk = create_group(loc, (*)(k...))
            kk["therapy"] = k
            kk["rebounds"] = outdict[k]
        end
    end
end

##
out1 = rebound_dict(ablist, 1 ) 
make_refinements!(out1;  quartile_estimator_samples = 20)

##
out2 = rebound_dict(ablist, 2 ) 
make_refinements!(out2;  quartile_estimator_samples = 20)

##
out3 = rebound_dict(ablist, 3 ) 
make_refinements!(out3;  quartile_estimator_samples = 20)

## Write output to h5 file

# Create file if it is not created
h5open(datapath*"/therapy_ranking.h5","w")
write_rebound_h5(out1)
write_rebound_h5(out2)
write_rebound_h5(out3)