*This repo is associated with https://github.com/CortXplorer/Dynamic_CSD*

# Dynamic CSD analysis steps in julia 

### AvrecPeakPipeline.jl

This pipeline takes extracted features from avrec and layer trace data in the form of csv (output from matlab) and produces graphs and plots for all t tests and cohen's d tests as well as a few extra plots to explore the data

Uses **AvrecPeakPlots.jl**, **AvrecPeakStats.jl**, and **CohensProg.jl** from /functions

Run ***RAnova.R*** for full anovas to complement these technically post-hoc tests. Bonferroni correction to posthoc t tests is n=14 p > 0.0038

### SpectralPermPipeline.jl

This pipeline takes .mat files from principle matlab CSD analysis and uses the single trial CSD data for broad spectral analysis and permutation testing (see Deane et al. 2020). It is a work in progress but doesn't look promising to continue.

Uses **callGroup.jl**, **CWTfunc.jl**, **Permfunctions.jl**, and **RunPermuations.jl**

### ThetaGammaCFC.jl

This script is a work in progress to run phase-amplitude cross frequency coupling analysis on theta vs low and high gamma bands

Uses **CFCfunc.jl**

### CSD_Pipeline.jl

This is a work in progress to start removing the need for Matlab from the conversion out of plexon rather than after initial csd analysis and avrec feature extraction. Currenlty to the point of generating a very small version of the .mat output in .jld2 from the raw converted files

Uses **csdStruct.jl**, **Dynamic_CSD.jl**, **SingleTrialCSD.jl**, **get_csd.jl**, **sink_dura.jl**, **functions.jl**, **callGroup.jl**