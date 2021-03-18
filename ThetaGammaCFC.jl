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
using FFTW, ProgressMeter, Random, PyCall # needed for just the getCFCcsv function
signal = pyimport("scipy.signal")         # python for the filtering
include(joinpath(func,"CFCfunc.jl"))      # contains supporting functions for CFC
include(joinpath(func,"getCFCcsv.jl"))    # contains main function for calculating CFC
sr = 1000                                 # sampling rate
NQ = Int(sr/2)                            # Nyquest frequency

getCFCcsv(home,data,Groups,condList,stimfrq,layers,sr,NQ)
#---------------------------------------------------------------------------------------

csvTab = CSV.File(joinpath(data,"Spectral\\CFCtable.csv")) |> DataFrame

# seperate by group
KIC = csvTab[csvTab[:,:Group] .== "KIC",:]
KIT = csvTab[csvTab[:,:Group] .== "KIT",:]
KIV = csvTab[csvTab[:,:Group] .== "KIV",:]

# performing this very roughly now to get a sense for it
KICpream   = KIC[KIC[:,:Condition] .== "preAM_1",:]
KICpream   = filter(row -> ! isnan(row.ObsDist_lowgam), KICpream) # find out where nans coming from!!
CprAM_mean = mean(KICpream.ObsDist_lowgam)

KITpream   = KIT[KIT[:,:Condition] .== "preAM_1",:]
KITpream   = filter(row -> ! isnan(row.ObsDist_lowgam), KITpream) # find out where nans coming from!!
TprAM_mean = mean(KITpream.ObsDist_lowgam)

KIVpream   = KIV[KIV[:,:Condition] .== "preAM_1",:]
KIVpream   = filter(row -> ! isnan(row.ObsDist_lowgam), KIVpream) # find out where nans coming from!!
VprAM_mean = mean(KIVpream.ObsDist_lowgam)

####

KICam   = KIC[KIC[:,:Condition] .== "AM_1",:]
KICam   = filter(row -> ! isnan(row.ObsDist_lowgam), KICam) # find out where nans coming from!!
CAM_mean = mean(KICam.ObsDist_lowgam)

KITam   = KIT[KIT[:,:Condition] .== "AM_1",:]
KITam   = filter(row -> ! isnan(row.ObsDist_lowgam), KITam) # find out where nans coming from!!
TAM_mean = mean(KITam.ObsDist_lowgam)

KIVam   = KIV[KIV[:,:Condition] .== "AM_1",:]
KIVam   = filter(row -> ! isnan(row.ObsDist_lowgam), KIVam) # find out where nans coming from!!
VAM_mean = mean(KIVam.ObsDist_lowgam)

####

KICprecl   = KIC[KIC[:,:Condition] .== "preCL_1",:];
KICprecl   = filter(row -> ! isnan(row.ObsDist_lowgam), KICprecl); # find out where nans coming from!!
CprCL_mean = mean(KICprecl.ObsDist_lowgam);

KITprecl   = KIT[KIT[:,:Condition] .== "preCL_1",:];
KITprecl   = filter(row -> ! isnan(row.ObsDist_lowgam), KITprecl); # find out where nans coming from!!
TprCL_mean = mean(KITprecl.ObsDist_lowgam);

KIVprecl   = KIV[KIV[:,:Condition] .== "preCL_1",:];
KIVprecl   = filter(row -> ! isnan(row.ObsDist_lowgam), KIVprecl); # find out where nans coming from!!
VprCL_mean = mean(KIVprecl.ObsDist_lowgam);

####

KICcl   = KIC[KIC[:,:Condition] .== "CL_1",:];
KICcl   = filter(row -> ! isnan(row.ObsDist_lowgam), KICcl); # find out where nans coming from!!
CCL_mean = mean(KICcl.ObsDist_lowgam);

KITcl   = KIT[KIT[:,:Condition] .== "CL_1",:];
KITcl   = filter(row -> ! isnan(row.ObsDist_lowgam), KITcl); # find out where nans coming from!!
TCL_mean = mean(KITcl.ObsDist_lowgam);

KIVcl   = KIV[KIV[:,:Condition] .== "CL_1",:];
KIVcl   = filter(row -> ! isnan(row.ObsDist_lowgam), KIVcl); # find out where nans coming from!!
VCL_mean = mean(KIVcl.ObsDist_lowgam);

####

println("Pre-laser AM: treated: ", TprAM_mean, ", naive and viral control: ", CprAM_mean, " and ", VprAM_mean)

println("Post-laser AM: treated: ", TAM_mean, ", naive and viral control: ", CAM_mean, " and ", VAM_mean)

println("Pre-laser CL: treated is ", TprCL_mean, ", naive and viral control: ", CprCL_mean, " and ", VprCL_mean)

println("Post-laser CL: treated: ", TCL_mean, ", naive and viral control: ", CCL_mean, " and ", VCL_mean)

using HypothesisTests, EffectSizes
TvCpream = pvalue(EqualVarianceTTest(KITpream.ObsDist_lowgam, KICpream.ObsDist_lowgam))
TvVpream = pvalue(EqualVarianceTTest(KITpream.ObsDist_lowgam, KIVpream.ObsDist_lowgam))
CvVpream = pvalue(EqualVarianceTTest(KICpream.ObsDist_lowgam, KIVpream.ObsDist_lowgam))

TvCam = pvalue(EqualVarianceTTest(KITam.ObsDist_lowgam, KICam.ObsDist_lowgam))
TvVam = pvalue(EqualVarianceTTest(KITam.ObsDist_lowgam, KIVam.ObsDist_lowgam))
CvVam = pvalue(EqualVarianceTTest(KICam.ObsDist_lowgam, KIVam.ObsDist_lowgam))

TvCprecl = pvalue(EqualVarianceTTest(KITprecl.ObsDist_lowgam, KICprecl.ObsDist_lowgam))
TvVprecl = pvalue(EqualVarianceTTest(KITprecl.ObsDist_lowgam, KIVprecl.ObsDist_lowgam))
CvVprecl = pvalue(EqualVarianceTTest(KICprecl.ObsDist_lowgam, KIVprecl.ObsDist_lowgam))

TvCcl = pvalue(EqualVarianceTTest(KITcl.ObsDist_lowgam, KICcl.ObsDist_lowgam))
TvVcl = pvalue(EqualVarianceTTest(KITcl.ObsDist_lowgam, KIVcl.ObsDist_lowgam))
CvVcl = pvalue(EqualVarianceTTest(KICcl.ObsDist_lowgam, KIVcl.ObsDist_lowgam))

TvCpream = effectsize(CohenD(KITpream.ObsDist_lowgam, KICpream.ObsDist_lowgam))
TvVpream = effectsize(CohenD(KITpream.ObsDist_lowgam, KIVpream.ObsDist_lowgam))
CvVpream = effectsize(CohenD(KICpream.ObsDist_lowgam, KIVpream.ObsDist_lowgam))

TvCam = effectsize(CohenD(KITam.ObsDist_lowgam, KICam.ObsDist_lowgam))
TvVam = effectsize(CohenD(KITam.ObsDist_lowgam, KIVam.ObsDist_lowgam))
CvVam = effectsize(CohenD(KICam.ObsDist_lowgam, KIVam.ObsDist_lowgam))

TvCprecl = effectsize(CohenD(KITprecl.ObsDist_lowgam, KICprecl.ObsDist_lowgam))
TvVprecl = effectsize(CohenD(KITprecl.ObsDist_lowgam, KIVprecl.ObsDist_lowgam))
CvVprecl = effectsize(CohenD(KICprecl.ObsDist_lowgam, KIVprecl.ObsDist_lowgam))

TvCcl = effectsize(CohenD(KITcl.ObsDist_lowgam, KICcl.ObsDist_lowgam))
TvVcl = effectsize(CohenD(KITcl.ObsDist_lowgam, KIVcl.ObsDist_lowgam))
CvVcl = effectsize(CohenD(KICcl.ObsDist_lowgam, KIVcl.ObsDist_lowgam))

## all pvalues significant, all cohen's d very small
## and there you have it. There doesn't seem to be a difference under anesthesia 