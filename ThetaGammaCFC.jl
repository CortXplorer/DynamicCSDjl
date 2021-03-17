# Theta Gamma Cross Frequency Coupling (CFC)

# Note: there are spots in the code specifically for extracting figures of the individual steps. These are currently silenced for the full loops but are a good sanity/understanding check 

# Code ported from: https://mark-kramer.github.io/Case-Studies-Python/07.html

using MAT, Statistics, DSP
using DataFrames, CSV

# set up directories
home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")
include(joinpath(group,"callGroup.jl")) # contains animal data

Groups   = ["KIC" "KIT" "KIV"]
condList = ["preCL_1","CL_1","preAM_1","AM_1"]
stimfrq  = ["2Hz" "5Hz" "10Hz" "20Hz" "40Hz"]
layers   = ["II" "IV" "V" "VI"]

# this calculates CFC at a single trial level and produces the CFCtable.csv with which further testing can be done. Current run-time is 8+ hours so ONLY run if a new table is needed. (we can figure out how to get this runtime down if needed)
# ______________________________________________________________________________________
# using FFTW, ProgressMeter, Random, PyCall # needed for just the getCFCcsv function
# signal = pyimport("scipy.signal")         # python for the filtering
# include(joinpath(func,"CFCfunc.jl"))      # contains supporting functions for CFC
# include(joinpath(func,"getCFCcsv.jl"))    # contains main function for calculating CFC
# sr = 1000                                 # sampling rate
# NQ = Int(sr/2)                            # Nyquest frequency
#
# getCFCcsv(home,data,signal,Groups,condList,stimfrq,layers,sr,NQ)
#---------------------------------------------------------------------------------------

csvTab = CSV.read(joinpath(data,"Spectral\\CFCtable.csv"))
