#!/usr/bin/env julia
# this script is thought to be designed to counduct SNAQ! on CFs
# which were obtaines from SNPs by SNPs2CF.
# input: 	- csv output of SNPs2CF
# 		- hypothetical tree as starting topology

#using Pkg # to use functions that manage packages
#Pkg.add("PhyloNetworks") # to download & install package PhyloNetworks
#Pkg.add("PhyloPlots")
#Pkg.add("RCall")      # packaage to call R from within julia
#Pkg.add("CSV")        # to read from / write to text files, e.g. csv files
#Pkg.add("DataFrames") # to create & manipulate data frames
#Pkg.add("StatsModels")# for regression formulas
using PhyloNetworks   # may take some time: pre-compiles functions in that package
# using PhyloPlots
# using Distributed

print("we are running")

snpCF = readTableCF(ARGS[1]) # input our CF table
tre = readTopology(ARGS[2]) # input a starting topology

# print(tre)

n0_out = ARGS[1][1:end-4] * "_net0"
# n1_out = ARGS[1][1:end-4] * "_net1"
# n2_out = ARGS[1][1:end-4] * "_net2"
# plot_out = ARGS[1][1:end-4] * "_plot.pdf"

# print(n0_out)

rseed = rand((100:999))
net0 = snaq!(tre,  snpCF, hmax=0, filename=n0_out, seed=rseed)
# rseed = rand((100:999))
# net1 = snaq!(net0, snpCF, hmax=1, filename=n1_out, seed=rseed)

# addprocs(2) # tells Julia to add 2 worker cores
# @everywhere using PhyloNetworks # tells these cores to load the PhyloNetworks package for themselves
# net2 = snaq!(net1, snpCF, hmax=2, filename=n2_out, seed=789, runs=5)

# plot(net0, :R);

# R"pdf"(plot_out, width=3, height=3);
# plot(net0, :R);
# R"dev.off()";


