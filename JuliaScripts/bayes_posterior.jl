include(joinpath(@__DIR__,"file_access.jl")) # run the 


##

#=
In this section we build the results of the analysis of the sites.
Single threaded and costly, takes about a day to run all the sites.
=#

filepath = datapath*"/snpanalysis.h5"
fid = h5open(filepath, "r+")
dictovec(a::Dict) = [a["$i"] for i in 1:length(a)] #auxiliary function for turning dictionaries into Vector
for ab in fid
    global bayesout = []
    for site in ab
        mut = reverse(read(site["mut"]))
        # A vector of vectors of the LikelihoodSample objects (structs with the random samples of the weight functions)
        ls = [ [
            LikelihoodSample(mut..., read(day["theta"]), dictovec(read(day["kwt_mut"]))... ; n_samp = 10000) 
            for day in patient if (sum(dictovec(read(day["kwt_mut"]))) > 0) & (read(day["theta"]) > 0)]
        for patient in site["snpdata"] ]
            for avg in [0, 1]
                posterior = baysian_pop_rs(ls; samples = 2*10^3, burn_in = 10^3, avg = avg)
                site["bayes_avg$(round(avg,digits=1))"] = posterior
                println("$avg") # print to see the progression through the sites
            end
    end
end
close(fid)

## Clade b limited version
#=

filepath = datapath*"/snpanalysis_cladeb.h5"
fid = h5open(filepath, "r+")
dictovec(a::Dict) = [a["$i"] for i in 1:length(a)] #auxiliary function for turning dictionaries into Vector
for ab in fid
    global bayesout = []
    for site in ab
        mut = reverse(read(site["mut"]))
        # A vector of vectors of the LikelihoodSample objects (structs with the random samples of the weight functions)
        ls = [ [
            LikelihoodSample(mut..., read(day["theta"]), dictovec(read(day["kwt_mut"]))... ; n_samp = 10000) 
            for day in patient if (sum(dictovec(read(day["kwt_mut"]))) > 0) & (read(day["theta"]) > 0)]
        for patient in site["snpdata"] ]
            for avg in [0, 1]
                posterior = baysian_pop_rs(ls; samples = 2*10^3, burn_in = 10^3, avg = avg)
                site["bayes_avg$(round(avg,digits=1))"] = posterior
                println("$avg") # print to see the progression through the sites
            end
    end
end
close(fid)

##


#$round(mutwt,digits = 2) & $round(mut[2],digits = 2) &  $round(r,digits=1)  $round(rσ, digits = 1) 

