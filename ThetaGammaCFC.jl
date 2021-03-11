# Theta Gamma Cross Frequency Coupling (CFC)

# https://mark-kramer.github.io/Case-Studies-Python/07.html

using MAT, Statistics
using DSP, FFTW
# python for the filtering:
using PyCall
signal = pyimport("scipy.signal")

# set up directories
home    = @__DIR__
data    = joinpath(home,"Data")
func    = joinpath(home,"functions")
include(joinpath(func,"CFCfunc.jl"))

# loop through groups
# KIT, KIC, KIV
# loop through animals
curAn   = "KIC02"
anDat   = matread(joinpath(data,(curAn * "_Data.mat")));
anCSD  = anDat["Data"]["SglTrl_CSD"]; # all single trial CSD data for full recording day
anCon  = anDat["Data"]["Condition"];  # all conditions for full recording day


# loop through relevant conditions
# let's start with the preCL
ConIdx = (LinearIndices(anCon))[findall(x -> x == "preCL_1", anCon)] #findall gives cartesian coordinates so we have to further translate to index position
curCSD = anCSD[ConIdx][1]; # pull out the CSD at that position
# loop through frequencies 
# loop through trials
thistrial = 6
curCSDst = curCSD[2][:,:,thistrial]  # pull out the 5hz signal for now in the FIRST single trial

# loop through layers (may flip order to layers then trials)
# to run through layers, check CWT_Loop from CWTfunc.jl, for now we select manually the middle channels of layer IV. That's our final signal to process:
layIV = curCSDst[7,:]

# Filter Step!
sr  = 1000          # sampling rate
NQ  = Int(sr/2)     # Nyquest frequency

# specplot = visSpectrum(layIV, sr) # to visualize the SPECTRUM

# set a passband [4-7] Hz 
Wn  = [4:7...] # theta band
n   = 100      # filter order (determines accuracy)
b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming") #PyCall
Vlo = signal.filtfilt(b, 1, layIV) #PyCall

# set a passband [31-60] Hz 
Wn  = [31:60...] #low gamma band (in KIC02 preCL_1 I'm seeing a spike at 50 Hz, I will have to find a way to verify that this is the target of interest accross measurements/animals)
b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming")
Vlg = signal.filtfilt(b, 1, layIV)

Wn  = [61:100...] #low gamma band (in KIC02 preCL_1 I'm seeing a spike at 50 Hz, I will have to find a way to verify that this is the target of interest accross measurements/animals)
b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming")
Vhg = signal.filtfilt(b, 1, layIV)

## visualize the filters
plot(
    [layIV Vlo Vlg Vhg], 
    label = ["Signal" "Theta" "Low Gamma" "High Gamma"],
    ylabel = "Amplitude [mV/mm²]",
    xlabel = "Time [ms]",
    title = "Spectrum",
)

phith = angle.(signal.hilbert(Vlo)) # Compute phase of low-freq signal
amplg = abs.(signal.hilbert(Vlg))   # Compute amplitude of high-freq signal
amphg = abs.(signal.hilbert(Vhg))   # Compute amplitude of high-freq signal

phasebins = [-pi:0.1:pi...] # To compute CFC, define phase bins,
amplgmean   = zeros(1,length(phasebins)-1) # preallocate 
amphgmean   = zeros(1,length(phasebins)-1) # preallocate 
phithmean   = zeros(1,length(phasebins)-1) # preallocate

for iBin = 1:length(phasebins)-1
    phLow = phasebins[iBin]   # lower limit
    phHi  = phasebins[iBin+1] # upper limit
    ind   = findall(phLow .<= phith .< phHi) # find phases falling in bin
    amplgmean[iBin] = mean(amplg[ind])         # compute mean amp
    amphgmean[iBin] = mean(amphg[ind])         # compute mean amp
    phithmean[iBin] = mean([phLow, phHi])    # save center phase 
end

### Determine if the Phase and Amp are related - of a variety of methods, they recommend the phase-amplitude plot:
plot([phithmean' phithmean'],
    [amplgmean' amphgmean'],
    linewidth = 4,
    linecolor = [:orange :red],
    title  = "Cross Frequency Coupling",
    xlabel = "Theta phase",
    ylabel = "Gamma amplitude",
    label = ["low gamma vs theta" "high gamma vs theta"],
    grid   = false,
)

hlg = maximum(amplgmean) - minimum(amplgmean) # diff between min and max
hhg = maximum(amphgmean) - minimum(amphgmean)


### Assess h's significance by creating a surrogate phase-amp vector by resampling without replacement the amplitude time series (explanation as to why in link above)

nsurrogate = 1000 # number of test surrogates 
hSlg  = zeros(nsurrogate)
for iSur = 1:nsurrogate 

    # n_surrogates = 1000;                    #Define no. of surrogates.
    # hS = zeros(n_surrogates)                #Vector to hold h results.
    # for ns in range(n_surrogates):          #For each surrogate,
    #     ampS = amp[randint(0,N,N)]          #Resample amplitude,
    #     p_bins = arange(-pi, pi, 0.1)       #Define the phase bins
    #     a_mean = zeros(size(p_bins)-1)      #Vector for average amps.
    #     p_mean = zeros(size(p_bins)-1)      #Vector for phase bins.
    #     for k in range(size(p_bins)-1):
    #         pL = p_bins[k]                  #... lower phase limit,
    #         pR = p_bins[k+1]                #... upper phase limit.
    #         indices=(phi>=pL) & (phi<pR)    #Find phases falling in bin,
    #         a_mean[k] = mean(ampS[indices]) #... compute mean amplitude,
    #         p_mean[k] = mean([pL, pR])      #... save center phase.
    #     hS[ns] = max(a_mean)-min(a_mean)    # Store surrogate h.

    # counts, _, _ = hist(hS, label='surrogates')               # Plot the histogram of hS, and save the bin counts.
    # vlines(h, 0, max(counts), colors='red', label='h', lw=2)  # Plot the observed h,
    # legend();
### - To compute a p-value, we determine the proportion of surrogate h values greater than the observed h value
    # p = sum([s > h for s in hS]) / len(hS)
    # print(p)