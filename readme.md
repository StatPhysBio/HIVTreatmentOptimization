The analysis pipelne assumes the stability of the file locations.
To recover the results:
1.) (Optional) Run the Mathematica files ending in "_analysis.nb", to (re)generate the HDF5 summary of
	the genetic and viremic results and estimate mut target size.
Make sure to run "snp_analysis.nb" first, since this file sets the site information
used in the genetic analysis.
In case you don't have access to Mathematica, 
	the .h5 files have been provided so you can test the Julia code 
	which follows (Julia is free).
Usually I run these sequentially in VSCode but it should be possible 
	to run them by calling:
`
$ julia file.jl
`
  
`bayes_posterior.jl` to generate posterior on the fiteness cost
	function. This is stored in the same file snpanalysis.h5
	
`fit_reservoir.jl` to perform the fitting 
	of the reservoir. You should be able to recover ξ = 2.1.
	
`simulate_therapy.jl` to simulate the result of viral rebounds and 
	compare optimal treatment.

`discrepancy_stats.jl` to generate the statistics used in the hypothesis testing

`linkage_simulations.jl` to test the inference procedure 
	(calls the Tomoko.jl package which runs longer genomes
	with symmetric mutation rates)

Tested with Julia 1.61.
 
