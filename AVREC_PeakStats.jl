using Statistics, DSP
using Plots, Colors
using CSV
using DataFrames

# Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
PeakData = CSV.File("AVRECPeakData.csv") |> DataFrame