# Theta Gamma Cross Frequency Coupling (CFC)

# Note: there are spots in the code specifically for extracting figures of the individual steps. These are currently silenced for the full loops but are a good sanity/understanding check 

# Code ported from: https://mark-kramer.github.io/Case-Studies-Python/07.html

using MAT, Statistics, Random
using DSP, FFTW
using DataFrames, CSV
# python for the filtering:
using PyCall
signal = pyimport("scipy.signal")

# set up directories
home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
group   = joinpath(home,"groups")
include(joinpath(func,"CFCfunc.jl")) # contains supporting functions for CFC
include(joinpath(group,"callGroup.jl")) # contains animal data

Groups   = ["KIC" "KIT" "KIV"]
condList = ["preCL_1","CL_1","preAM_1","AM_1"]
stimfrq  = ["2Hz" "5Hz" "10Hz" "20Hz" "40Hz"]
layers   = ["II" "IV" "V" "VI"]
sr       = 1000          # sampling rate
NQ       = Int(sr/2)     # Nyquest frequency

# initialize DataFrame table
cfcTab = DataFrame(Group = String[], Animal = String[], Condition = String[], StimFrq = String[], Layer = String[], h_lowgam = Float64[], Smean_lowgam = Float64[], Sstd_lowgam = Float64[], p_lowgam = Float64[], ObsDist_lowgam = Float64[], h_higam = Float64[], Smean_higam = Float64[], Sstd_higam = Float64[], p_higam = Float64[], ObsDist_higam = Float64[])

takepic = 0 # only use this for specific case checking

for iGr = 1:length(Groups) # loop through groups

    animalList,_,LIIList,LIVList,LVList,LVIList,_ = callGroup(Groups[iGr]); # extract animal data of group

    for iAn = 1:length(animalList) # loop through animals

        curAn   = animalList[iAn]
        anDat   = matread(joinpath(data,(curAn * "_Data.mat")));
        anCSD  = anDat["Data"]["SglTrl_CSD"]; # all single trial CSD data 
        anCon  = anDat["Data"]["Condition"];  # all conditions for full recording day

        for iCn = 1:length(condList) # loop through relevant conditions

            ConIdx = findall(x -> x == condList[iCn], anCon) 
            curCSD = anCSD[ConIdx][1]; # pull out the CSD at that condition

            for iFr = 1:length(stimfrq) # loop through frequencies 
            
                for iTr = 1:size(curCSD[iFr],3) # loop through trials
                    
                    curCSDst = curCSD[iFr][:,:,iTr]  # pull out the signal

                    for iLa = 1:length(layers) # loop through layers 
                    # to run through layers, check CWT_Loop from CWTfunc.jl, for now we select manually the middle channels of layer IV. That's our final signal to process:

                        if iLa == 1 # we can write this to check correct layers later
                            curChan = LIIList[iAn]
                        elseif iLa == 2
                            curChan = LIVList[iAn]
                        elseif iLa == 3
                            curChan = LVList[iAn]
                        elseif iLa == 4
                            curChan = LVIList[iAn]
                        end

                        centerChan = Int(ceil(mean(curChan))) # take only the central layer channel

                        layCSDst = curCSDst[centerChan,:] # now we have our final raw signal to start with!

                        # Filter Step!
                        Vlo, Vlg, Vhg = thetagamma_filter(layCSDst, sr, NQ, signal, takepic)

                        phith = angle.(signal.hilbert(Vlo)) # Compute phase of theta
                        amplg = abs.(signal.hilbert(Vlg))   # Compute amp of low gamma
                        amphg = abs.(signal.hilbert(Vhg))   # Compute amp of high gamma

                        hlg, hhg = CFCget_h(phith,amplg,amphg,1) # h gives the max amp - min amp after creating the curve to couple theta phase to gamma amp

                        ### Assess h's significance by creating a surrogate phase-amp vector by resampling without replacement the amplitude time series (explanation as to why in link above)

                        nsurrogate = 1000 # number of test surrogates 
                        hSlg  = zeros(nsurrogate) # very like a parametric bootstrap
                        hShg  = zeros(nsurrogate) # also very like permuation testing

                        for iSur = 1:nsurrogate 

                            amplgS = shuffle(amplg) # shuffle or resample
                            amphgS = shuffle(amphg)

                            hSlg[iSur], hShg[iSur] = CFCget_h(phith,amplgS,amphgS,0)
                            
                        end

                        if takepic == 1
                            distplot = visSurrdist(hSlg, hShg, hlg, hhg)
                            savefig("distSurr_vsObs.png")
                        end

                        ### - To compute a p-value, we determine the proportion of surrogate h values greater than the observed h value

                        pval_lg = length(findall(hSlg .> hlg)) / length(hSlg) # p value
                        pval_hg = length(findall(hShg .> hhg)) / length(hShg)

                        surmean_lg, surstd_lg = mean(hSlg), std(hSlg) # mean/std of sur dist
                        distObs_lg = hlg - surmean_lg / surstd_lg     # std away from mean 

                        surmean_hg, surstd_hg = mean(hShg), std(hShg)
                        distObs_hg = hhg - surmean_hg / surstd_hg

                        push!(cfcTab, [curAn[1:3] curAn condList[iCn] stimfrq[iFr] layers[iLa] hlg surmean_lg surstd_lg pval_lg distObs_lg hhg surmean_hg surstd_hg pval_hg distObs_hg])

                        # averaging the std distance of observed from the surrogate mean between trials and then animals will give a magnitude of significant difference from random and may provide more nuanced insights into differences between groups 
                    end # layers
                end # trials
            end # frequency
        end # conditions
    end # animals
end # groups

CSV.write(joinpath(data,"Spectral\\CFCtable.csv"), cfcTab) # write out dataframe to CSV