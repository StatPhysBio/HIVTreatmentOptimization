The analysis pipelne assumes the stability of the file locations.
To recover the results:
1.) (Optional) Run the Mathematica files 
	in GeneticViremicAnalysis to (re)generate the HDF5 summary of
	the genetic and viremic results and estimate mut target size.
Make sure to run "snp_analysis.nb" first, since this file sets the site information
used in the genetic analysis.
In case you don't have access to Mathematica, 
	the .h5 files have been provided so you can test the Julia code 
	which follows (Julia is free).
2.) Run Julia scripts to do the heavy lifting. I run these sequentially in VSCode but it should be possible 
	to run them from the command line.
	These call packages EscapeSimulator, and the auxiliary functions defined in `file_access.jl` that make it easy to simulate using fitted paramters.
  
`bayes_posterior.jl` to generate posterior on the fiteness cost
	function. This is stored in the same file snpanalysis.h5
	
`fit_reservoir.jl` to perform the fitting 
	of the reservoir. You should be able to recover ξ = 2.1.
	Also generates the statistics used in the hypothesis testing
	
`simulate_therapy.jl` to simulate the 
	result of viral rebounds in trials

`therapy_ranking.jl` to compare and rank therapies with arbitrary 			number of combinations of antibodies.

`linkage_simulations.jl` to test the inference procedure 
	(calls the Tomoko.jl package which runs longer genomes
	with symmetric mutation rates)

Tested with Julia 1.61.
 3.) Finally we include some (Mathematica-based) plotting functions in the Plotting folder. This includes a file of helper functions `PlotFunctions.nb`

and a file for figure generation.

 `PlotGeneration.nb`

Please feel free to file an issues if something is confusing or needs explanation.