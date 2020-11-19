# compute scalograms 
using MAT
using JLD

home    = @__DIR__
path    = joinpath(home,"Data")
anDat   = matread(joinpath(path,"KIC02_Data.mat"))
jldDat  = load(joinpath(path,"KIT02_Data.jld"))

var = 3
name = joinpath(path,"vari3.jld")
save(name, "t", t, "arr", z)


# generate spectral plots (osci bands sorted by group and then frequency)

# spectral power permutations
# between group
# within group

# spectral phase coherence permutations
# between group
# within group