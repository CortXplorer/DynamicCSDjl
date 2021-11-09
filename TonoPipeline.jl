#____________________________________________________________________________________
# Feature extraction of best frequency from all tonotopies listed 
#____________________________________________________________________________________
using MAT, StatsBase, DSP
using DataFrames, CSV
using Peaks

home    = @__DIR__
groups  = joinpath(home,"groups")
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")

include(joinpath(groups,"callGroup.jl"))
include(joinpath(func,"TonotopyAnalysis.jl"))
include(joinpath(func,"nanmean.jl"))
include(joinpath(func,"ifindpeaks.jl"))

# function takes all tonotopies in current data (specific call by measurement names inside function) and performs feature extraction to compare with click and AM features 
# input:    .../Data/*_Data.mat; .../groups/callGroup.jl
# output:   .../Data/Tonotopy/TonoRMS.csv
TonotopyAnalysis(data)  

#____________________________________________________________________________________
# Progression of strength of response (RMS) over course of experiment
#____________________________________________________________________________________

using StatsPlots, DataFrames, CSV

home    = @__DIR__
data    = joinpath(home,"Data")
name    = joinpath(data,"Tonotopy","TonoRMS") * ".csv"
TonoRMS = CSV.File(name) |> DataFrame

TonoProgressionAVREC = @df TonoRMS groupedboxplot(:Measurement,:AVRECrms, group = :Group, markerstrokewidth=0, legend = true, xrotation = 45)

TonoRMSIV = filter(row -> ! isnan(row.LIVrms), TonoRMS)

TonoProgressionLIV = @df TonoRMSIV groupedboxplot(:Measurement,:LIVrms, group = :Group, markerstrokewidth=0, legend = true, xrotation = 45)

# STATISTICS, difference between groups for each tonotopy measurement 

#____________________________________________________________________________________
# Peak amp vs lat for ALL measurement types
#____________________________________________________________________________________

using StatsPlots, DataFrames, CSV
using LsqFit

home    = @__DIR__
data    = joinpath(home,"Data")
name    = joinpath(data,"Tonotopy","TonoRMS") * ".csv"
TonoRMS = CSV.File(name) |> DataFrame
name    = joinpath(data,"AVRECPeakAMST") * ".csv"
AvrecAM = CSV.File(name) |> DataFrame
name    = joinpath(data,"AVRECPeakCLST") * ".csv"
AvrecCL = CSV.File(name) |> DataFrame

# just first click
AvrecCL = AvrecCL[AvrecCL[:,:OrderofClick] .== 1,:]
AvrecAM = AvrecAM[AvrecAM[:,:OrderofClick] .== 1,:]
# just 5 Hz
AvrecCL = AvrecCL[AvrecCL[:,:ClickFreq] .== 5,:]
AvrecAM = AvrecAM[AvrecAM[:,:ClickFreq] .== 5,:]

# combine dataframes, plot distributions? exp function? 
# AVREC
# CL data
TAllCL = AvrecCL[AvrecCL[:,:Group] .== "KIT",:]
TAllCL = TAllCL[TAllCL[:,:Layer] .== "All",:]
TAllCL = filter(row -> ! isnan(row.PeakAmp), TAllCL)
Tclamp = TAllCL[:,:PeakAmp]
Tcllat = TAllCL[:,:PeakLat]
# AM data
TAllAM = AvrecAM[AvrecAM[:,:Group] .== "KIT",:]
TAllAM = TAllAM[TAllAM[:,:Layer] .== "All",:]
TAllAM = filter(row -> ! isnan(row.PeakAmp), TAllAM)
Tamamp = TAllAM[:,:PeakAmp]
Tamlat = TAllAM[:,:PeakLat]
# Tono data
TAllTo = TonoRMS[TonoRMS[:,:Group] .== "KIT",:]
TAllTo = filter(row -> ! isnan(row.AVRECamp), TAllTo)
Ttoamp = TAllTo[:,:AVRECamp]
Ttolat = TAllTo[:,:AVREClat]

Tamp = cat(Tclamp,Tamamp,Ttoamp,dims=1)
Tlat = cat(Tcllat,Tamlat,Ttolat,dims=1)
# scatter(Tlat,Tamp,alpha=0.1)

xval = sort(unique(Tlat))
avg  = zeros(length(xval))
ste  = zeros(length(xval))
for ix = 1:length(xval)
    avg[ix] = mean(Tamp[Tlat.==xval[ix]])
    ste[ix] = sem(Tamp[Tlat.==xval[ix]])
end

plot!(xval,avg)

# CL data
CAllCL = AvrecCL[AvrecCL[:,:Group] .== "KIC",:]
CAllCL = CAllCL[CAllCL[:,:Layer] .== "All",:]
CAllCL = filter(row -> ! isnan(row.PeakAmp), CAllCL)
Cclamp = CAllCL[:,:PeakAmp]
Ccllat = CAllCL[:,:PeakLat]
# AM data
CAllAM = AvrecAM[AvrecAM[:,:Group] .== "KIC",:]
CAllAM = CAllAM[CAllAM[:,:Layer] .== "All",:]
CAllAM = filter(row -> ! isnan(row.PeakAmp), CAllAM)
Camamp = CAllAM[:,:PeakAmp]
Camlat = CAllAM[:,:PeakLat]
# Tono data
CAllTo = TonoRMS[TonoRMS[:,:Group] .== "KIC",:]
CAllTo = filter(row -> ! isnan(row.AVRECamp), CAllTo)
Ctoamp = CAllTo[:,:AVRECamp]
Ctolat = CAllTo[:,:AVREClat]

Camp = cat(Cclamp,Camamp,Ctoamp,dims=1)
Clat = cat(Ccllat,Camlat,Ctolat,dims=1)
# scatter!(Clat,Camp,alpha=0.1)

xval = sort(unique(Clat))
avg  = zeros(length(xval))
ste  = zeros(length(xval))
for ix = 1:length(xval)
    avg[ix] = mean(Camp[Clat.==xval[ix]])
    ste[ix] = sem(Camp[Clat.==xval[ix]])
end

plot!(xval,avg)

# CL data
VAllCL = AvrecCL[AvrecCL[:,:Group] .== "KIV",:]
VAllCL = VAllCL[VAllCL[:,:Layer] .== "All",:]
VAllCL = filter(row -> ! isnan(row.PeakAmp), VAllCL)
Vclamp = VAllCL[:,:PeakAmp]
Vcllat = VAllCL[:,:PeakLat]
# AM data
VAllAM = AvrecAM[AvrecAM[:,:Group] .== "KIV",:]
VAllAM = VAllAM[VAllAM[:,:Layer] .== "All",:]
VAllAM = filter(row -> ! isnan(row.PeakAmp), VAllAM)
Vamamp = VAllAM[:,:PeakAmp]
Vamlat = VAllAM[:,:PeakLat]
# Tono data
VAllTo = TonoRMS[TonoRMS[:,:Group] .== "KIV",:]
VAllTo = filter(row -> ! isnan(row.AVRECamp), VAllTo)
Vtoamp = VAllTo[:,:AVRECamp]
Vtolat = VAllTo[:,:AVREClat]

Vamp = cat(Vclamp,Vamamp,Vtoamp,dims=1)
Vlat = cat(Vcllat,Vamlat,Vtolat,dims=1)
# scatter!(Vlat,Vamp,alpha=0.1)

xval = sort(unique(Vlat))
avg  = zeros(length(xval))
ste  = zeros(length(xval))
for ix = 1:length(xval)
    avg[ix] = mean(Vamp[Vlat.==xval[ix]])
    ste[ix] = sem(Vamp[Vlat.==xval[ix]])
end

plot!(xval,avg)
 

model(t, p) = p[1] * exp.(-p[2] * t)
p0 = [0.5, 0.5]
fit = curve_fit(model, Tlat, Tamp, p0)
param = fit.param
cov = estimate_covar(fit)
lineT = model(Tlat,param)