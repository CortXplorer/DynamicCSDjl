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
stimfrq  = ["5Hz"] #"2Hz" "5Hz" "10Hz" "20Hz" "40Hz"
layers   = ["IV" "V"] # "II" "IV" "V" "VI"

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

KIC = KIC[KIC[:,:Layer] .== "VI",:]
KIT = KIT[KIT[:,:Layer] .== "VI",:]
KIV = KIV[KIV[:,:Layer] .== "VI",:]

# performing this very roughly now to get a sense for it
KICpream   = KIC[KIC[:,:Condition] .== "preAM_1",:]
KICpream   = filter(row -> ! isnan(row.ObsDist_higam), KICpream) # find out where nans coming from!!
CprAM_mean = mean(KICpream.ObsDist_higam)

KITpream   = KIT[KIT[:,:Condition] .== "preAM_1",:]
KITpream   = filter(row -> ! isnan(row.ObsDist_higam), KITpream) # find out where nans coming from!!
TprAM_mean = mean(KITpream.ObsDist_higam)

KIVpream   = KIV[KIV[:,:Condition] .== "preAM_1",:]
KIVpream   = filter(row -> ! isnan(row.ObsDist_higam), KIVpream) # find out where nans coming from!!
VprAM_mean = mean(KIVpream.ObsDist_higam)

####

KICam   = KIC[KIC[:,:Condition] .== "AM_1",:]
KICam   = filter(row -> ! isnan(row.ObsDist_higam), KICam) # find out where nans coming from!!
CAM_mean = mean(KICam.ObsDist_higam)

KITam   = KIT[KIT[:,:Condition] .== "AM_1",:]
KITam   = filter(row -> ! isnan(row.ObsDist_higam), KITam) # find out where nans coming from!!
TAM_mean = mean(KITam.ObsDist_higam)

KIVam   = KIV[KIV[:,:Condition] .== "AM_1",:]
KIVam   = filter(row -> ! isnan(row.ObsDist_higam), KIVam) # find out where nans coming from!!
VAM_mean = mean(KIVam.ObsDist_higam)

####

KICprecl   = KIC[KIC[:,:Condition] .== "preCL_1",:];
KICprecl   = filter(row -> ! isnan(row.ObsDist_higam), KICprecl); # find out where nans coming from!!
CprCL_mean = mean(KICprecl.ObsDist_higam);

KITprecl   = KIT[KIT[:,:Condition] .== "preCL_1",:];
KITprecl   = filter(row -> ! isnan(row.ObsDist_higam), KITprecl); # find out where nans coming from!!
TprCL_mean = mean(KITprecl.ObsDist_higam);

KIVprecl   = KIV[KIV[:,:Condition] .== "preCL_1",:];
KIVprecl   = filter(row -> ! isnan(row.ObsDist_higam), KIVprecl); # find out where nans coming from!!
VprCL_mean = mean(KIVprecl.ObsDist_higam);

####

KICcl   = KIC[KIC[:,:Condition] .== "CL_1",:];
KICcl   = filter(row -> ! isnan(row.ObsDist_higam), KICcl); # find out where nans coming from!!
CCL_mean = mean(KICcl.ObsDist_higam);

KITcl   = KIT[KIT[:,:Condition] .== "CL_1",:];
KITcl   = filter(row -> ! isnan(row.ObsDist_higam), KITcl); # find out where nans coming from!!
TCL_mean = mean(KITcl.ObsDist_higam);

KIVcl   = KIV[KIV[:,:Condition] .== "CL_1",:];
KIVcl   = filter(row -> ! isnan(row.ObsDist_higam), KIVcl); # find out where nans coming from!!
VCL_mean = mean(KIVcl.ObsDist_higam);

####

println("Pre-laser AM: treated: ", TprAM_mean, ", naive and viral control: ", CprAM_mean, " and ", VprAM_mean)

println("Post-laser AM: treated: ", TAM_mean, ", naive and viral control: ", CAM_mean, " and ", VAM_mean)

println("Pre-laser CL: treated is ", TprCL_mean, ", naive and viral control: ", CprCL_mean, " and ", VprCL_mean)

println("Post-laser CL: treated: ", TCL_mean, ", naive and viral control: ", CCL_mean, " and ", VCL_mean)

using HypothesisTests, EffectSizes
TvCpream = pvalue(EqualVarianceTTest(KITpream.ObsDist_higam, KICpream.ObsDist_higam))
TvVpream = pvalue(EqualVarianceTTest(KITpream.ObsDist_higam, KIVpream.ObsDist_higam))
CvVpream = pvalue(EqualVarianceTTest(KICpream.ObsDist_higam, KIVpream.ObsDist_higam))

TvCam = pvalue(EqualVarianceTTest(KITam.ObsDist_higam, KICam.ObsDist_higam))
TvVam = pvalue(EqualVarianceTTest(KITam.ObsDist_higam, KIVam.ObsDist_higam))
CvVam = pvalue(EqualVarianceTTest(KICam.ObsDist_higam, KIVam.ObsDist_higam))

TvCprecl = pvalue(EqualVarianceTTest(KITprecl.ObsDist_higam, KICprecl.ObsDist_higam))
TvVprecl = pvalue(EqualVarianceTTest(KITprecl.ObsDist_higam, KIVprecl.ObsDist_higam))
CvVprecl = pvalue(EqualVarianceTTest(KICprecl.ObsDist_higam, KIVprecl.ObsDist_higam))

TvCcl = pvalue(EqualVarianceTTest(KITcl.ObsDist_higam, KICcl.ObsDist_higam))
TvVcl = pvalue(EqualVarianceTTest(KITcl.ObsDist_higam, KIVcl.ObsDist_higam))
CvVcl = pvalue(EqualVarianceTTest(KICcl.ObsDist_higam, KIVcl.ObsDist_higam))

TvCpreamEF = effectsize(CohenD(KITpream.ObsDist_higam, KICpream.ObsDist_higam))
TvVpreamEF = effectsize(CohenD(KITpream.ObsDist_higam, KIVpream.ObsDist_higam))
CvVpreamEF = effectsize(CohenD(KICpream.ObsDist_higam, KIVpream.ObsDist_higam))

TvCamEF = effectsize(CohenD(KITam.ObsDist_higam, KICam.ObsDist_higam))
TvVamEF = effectsize(CohenD(KITam.ObsDist_higam, KIVam.ObsDist_higam))
CvVamEF = effectsize(CohenD(KICam.ObsDist_higam, KIVam.ObsDist_higam))

TvCpreclEF = effectsize(CohenD(KITprecl.ObsDist_higam, KICprecl.ObsDist_higam))
TvVpreclEF = effectsize(CohenD(KITprecl.ObsDist_higam, KIVprecl.ObsDist_higam))
CvVpreclEF = effectsize(CohenD(KICprecl.ObsDist_higam, KIVprecl.ObsDist_higam))

TvCclEF = effectsize(CohenD(KITcl.ObsDist_higam, KICcl.ObsDist_higam))
TvVclEF = effectsize(CohenD(KITcl.ObsDist_higam, KIVcl.ObsDist_higam))
CvVclEF = effectsize(CohenD(KICcl.ObsDist_higam, KIVcl.ObsDist_higam))

## all pvalues significant, all cohen's d very small or small

Tprevcl = pvalue(EqualVarianceTTest(KITprecl.ObsDist_higam, KITcl.ObsDist_higam))
TprevclEF = effectsize(CohenD(KITprecl.ObsDist_higam, KITcl.ObsDist_higam))

Tprevam = pvalue(EqualVarianceTTest(KITpream.ObsDist_higam, KITam.ObsDist_higam))
TprevamEF = effectsize(CohenD(KITpream.ObsDist_higam, KITam.ObsDist_higam))

## within group before and after laser for treated shows no difference
