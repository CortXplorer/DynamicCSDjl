# Theta Gamma Cross Frequency Coupling (CFC)

# https://mark-kramer.github.io/Case-Studies-Python/07.html

using MAT, Statistics
using DSP, FFTW

home    = @__DIR__
data    = joinpath(home,"Data")

curAn   = "KIC02"
anDat   = matread(joinpath(data,(curAn * "_Data.mat")));

# all single trial CSD data for full recording day:
anCSD  = anDat["Data"]["SglTrl_CSD"];
# all conditions for recording day with this animal:
anCon  = anDat["Data"]["Condition"];

# let's start with the preCLtono
ConIdx = (LinearIndices(anCon))[findall(x -> x == "preCL_1", anCon)] #findall gives cartesian coordinates so we have to further translate to index position
# pull out the CSD at that position
curCSD = anCSD[ConIdx][1];
# pull out the 5 hz signal for now as well as the FIRST single trial
curCSD = curCSD[2][:,:,1]

# to run through layers, check CWT_Loop from CWTfunc.jl, for now we select manually the middle 3 channels of layer IV. Then we average them and that's our final signal to process
layIV = mean(curCSD[6:8,:], dims =1)
# extract time, sampling interval, and Nyquist frequency
t   = [1:1:length(layIV)...] #theoretically to account for other sampling rates we can make (1000/sr) which would give us "1" and someone with a lower sampling rate a larger step size
dt  = t[2] - t[1] #sampling rate is 1000
fNQ = 1 / dt / 2
N   = length(layIV) # number of data points 
T   = length(layIV) # the total duration of data (T = N here)

# ## let's visualize the SPECTRUM
# x  = hanning(N) .* layIV[:] # multiply data by hanning taper
# xf = fft(x) # compute Fourier transform
# Sxx = (2*dt)^(2/T) * (xf.*conj(xf)) # compute the spectrum
# Sxx = real(Sxx) # ignore complex components
# df = 1/T # define frequency resolution 
# faxis = [0:df:fNQ+df...] # construct freq axis 
# plot(
#     faxis[1:139], 
#     10 * log10.(Sxx[1:139]),
#     title  = "Spectrum of CSD FFT",
#     ylabel = "Power [mV^2/Hz]",
#     xlabel = "Frequency [Hz]", 
# )

# python is better at this filtering stuff
using PyCall
signal = pyimport("scipy.signal")
# set a passband [4-7] Hz 
Wn  = [4:7...] # theta band
n   = 100
b   = signal.firwin(n, Wn, nyq=500, pass_zero=false, window="hamming")
Vlo = signal.filtfilt(b, 1, layIV)

# set a passband [31-60] Hz 
Wn  = [31:60...] #low gamma band (in KIC02 I'm seeing a spike at 50 Hz)
b   = signal.firwin(n, Wn, nyq=500, pass_zero=false, window="hamming")
Vhi = signal.filtfilt(b, 1, layIV)

## visualize the filters
plot(layIV', label = "signal")
plot!(Vlo', label = "theta")
plot!(Vhi', label = "low gamma")

    # phi = angle(signal.hilbert(Vlo))     # Compute phase of low-freq signal
    # amp = abs(signal.hilbert(Vhi))       # Compute amplitude of high-freq signal

    # p_bins = arange(-pi,pi,0.1)          # To compute CFC, define phase bins,
    # a_mean = zeros(size(p_bins)-1)       # ... variable to hold the amplitude,
    # p_mean = zeros(size(p_bins)-1)       # ... and variable to hold the phase.
    # for k in range(size(p_bins)-1):      # For each phase bin,
    #     pL = p_bins[k]                   #... get lower phase limit,
    #     pR = p_bins[k+1]                 #... get upper phase limit.
    #     indices=(phi>=pL) & (phi<pR)     #Find phases falling in this bin,
    #     a_mean[k] = mean(amp[indices])   #... compute mean amplitude,
    #     p_mean[k] = mean([pL, pR])       #... save center phase.
    # plot(p_mean, a_mean)                 #Plot the phase versus amplitude,
    # ylabel('High-frequency amplitude')   #... with axes labeled.
    # xlabel('Low-frequency phase')
    # title('CFC');


## Step by step: 
### Compute Spectrum of LFP data - analyze entire length of data and compute the spectrum with a Hanning taper
    # dt = t[1] - t[0]                # Define the sampling interval,
    # T = t[-1]                       # ... the duration of the data,
    # N = len(LFP)                    # ... and the no. of data points

    # x = hanning(N) * LFP            # Multiply data by a Hanning taper
    # xf = rfft(x - x.mean())         # Compute Fourier transform
    # Sxx = 2*dt**2/T * (xf*conj(xf)) # Compute the spectrum
    # Sxx = real(Sxx)                 # Ignore complex components

    # df = 1 / T                      # Define frequency resolution,
    # fNQ = 1 / dt / 2                # ... and Nyquist frequency. 

    # faxis = arange(0, fNQ + df, df) # Construct freq. axis
    # plot(faxis, 10 * log10(Sxx))    # Plot spectrum vs freq.
    # xlim([0, 200])                  # Set freq. range, 
    # ylim([-80, 0])                  # ... and decibel range
    # xlabel('Frequency [Hz]')        # Label the axes
    # ylabel('Power [mV$^2$/Hz]');

### Filter Data into High and Low frequency bands - should be two frequency bands chosen because they are the two highest peaks in the spectrum. The following is with the Finite Impulse Response (FIR) filter
    # from scipy import signal
    # Wn = [5,7];                         # Set the passband [5-7] Hz,
    # n = 100;                            # ... and filter order,
    #                                     # ... build the bandpass filter,
    # b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
    # Vlo = signal.filtfilt(b, 1, LFP);   # ... and apply it to the data.

    # Wn = [80, 120];                     # Set the passband [80-120] Hz,
    # n = 100;                            # ... and filter order,
    #                                     # ... build the bandpass filter,
    # b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
    # Vhi = signal.filtfilt(b, 1, LFP);   # ... and apply it to the data.
    # figure(figsize=(14, 4))         # Create a figure with a specific size.
    # plot(t, LFP)                    # Plot the original data vs time.
    # plot(t, Vlo)                    # Plot the low-frequency filtered data vs time.
    # plot(t, Vhi)                    # Plot the high-frequency filtered data vs time.
    # xlabel('Time [s]')
    # xlim([24, 26]);                 # Choose a 2 s interval to examine
    # ylim([-2, 2]);
    # legend(['LFP', 'Vlo', 'Vhi']);  # Add a legend.

### Extract amp and phase from filtered signals - phase from low-frequency signal and amp envelope of high-freq signal. They use the analytic signal approach to estimate instantaneous phase and amplitude envelope of the LFP, using a hilbert transform (lengthy description of why in link above)
    # phi = angle(signal.hilbert(Vlo))  # Compute phase of low-freq signal
    # amp = abs(signal.hilbert(Vhi))       # Compute amplitude of high-freq signal

### Determine if the Phase and Amp are related - of a variety of methods, they recommend the phase-amplitude plot
    # p_bins = arange(-pi, pi, 0.1)
    # a_mean = zeros(size(p_bins)-1)
    # p_mean = zeros(size(p_bins)-1)
    # for k in range(size(p_bins)-1):     #For each phase bin,
    #     pL = p_bins[k]                  #... lower phase limit,
    #     pR = p_bins[k+1]                #... upper phase limit.
    #     indices=(phi>=pL) & (phi<pR)    #Find phases falling in bin,
    #     a_mean[k] = mean(amp[indices])  #... compute mean amplitude,
    #     p_mean[k] = mean([pL, pR])      #... save center phase.
    # plot(p_mean, a_mean)                #Plot the phase versus amplitude,
    # ylabel('High-frequency amplitude')  #... with axes labeled.
    # xlabel('Low-frequency phase');
### - find the difference between the max and min of the average amp over phases
    # h = max(a_mean)-min(a_mean)
    # print(h)
### - assess h's significance by creating a surrogate phase-amp vector by resampling without replacement the amplitude time series (explanation as to why in link above)
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