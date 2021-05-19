function visSpectrum(layLFPst, sr = 1000)

    t   = [1:Int(1000/sr):length(layLFPst)...] #time accounting for sr
    dt  = t[2] - t[1]   # sampling interval
    fNQ = 1 / dt / 2    # Nyquest is half sampling rate (fraction: 0.5)
    N   = length(t)     # number of data points 
    T   = length(layLFPst) # the total duration of data (T = N here)

    x  = hanning(N) .* layLFPst[:] # multiply data by hanning taper
    xf = fft(x) # compute Fourier transform
    Sxx = (2*dt)^(2/T) * (xf.*conj(xf)) # compute the spectrum
    Sxx = real(Sxx) # ignore complex components
    df = 1/T # define frequency resolution 
    faxis = [0:df:fNQ+df...] # construct freq axis 
    specplot = plot(
        faxis[1:139], 
        10 * log10.(Sxx[1:139]),
        title  = "Spectrum of LFP FFT",
        ylabel = "Power [mV^2/Hz]",
        xlabel = "Frequency [Hz]", 
    )

    return specplot

end

function visSurrdist(hSlg, hShg, hlg, hhg)
    distplot = histogram(
        [hSlg hShg],
        label = ["Low gamma surrogate distribution" "High gamma surrogate distribution"],    
        fillcolor = [:red4 :blue4],
        title     = "Observed vs Surrogate Distribution",
        ylabel    = "Count",
        xlabel    = "Amplitude [mV/mm²]"
    ) # the results of the surrogate tests in a distribution

    vline!(
        [hlg hhg],
        label     = ["Low gamma observed" "High gamma observed"], 
        linecolor = [:indianred :royalblue2],
        linewidth = 6
    ) # with the observed values overlaid

    return distplot 

end

function thetagamma_filter(layLFPst, sr, NQ, signal, takepic=1)

    # set a passband [4-7] Hz 
    Wn  = [4:7...] # theta band
    n   = 100      # filter order (determines accuracy)
    b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming") #PyCall
    Vlo = signal.filtfilt(b, 1, layLFPst) #PyCall

    # set a passband [31-60] Hz 
    Wn  = [31:60...] #low gamma band (in KIC02 preCL_1 I'm seeing a spike at 50 Hz, I will have to find a way to verify that this is the target of interest accross measurements/animals)
    b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming")
    Vlg = signal.filtfilt(b, 1, layLFPst)

    Wn  = [61:100...] #high gamma band (I don't see a spike here but I figured we should still just check it for now)
    b   = signal.firwin(n, Wn, nyq=NQ, pass_zero=false, window="hamming")
    Vhg = signal.filtfilt(b, 1, layLFPst)

    if takepic == 1 # if we want mass pic production, will need to update naming labels
        ## visualize the spectrum
        specplot = visSpectrum(layLFPst, sr) # to visualize the SPECTRUM
        savefig("VisSpectrum.png") # yea please make a folder for this output, jesus
        # this is relevant in selecting the filter sizes and will need to be checked again later to verify our choices above

        ## visualize the filters
        plot(
            [layLFPst Vlo Vlg Vhg], 
            label = ["Signal" "Theta" "Low Gamma" "High Gamma"],
            ylabel = "Amplitude [mV/mm²]",
            xlabel = "Time [ms]",
            title = "Filters",
        )
        savefig("VisFilters.png") # also folder for this
    end



    return Vlo, Vlg, Vhg
end

function CFCget_h(phith,amplg,amphg,takepic=1)

    phasebins = [-pi:0.1:pi...] # To compute CFC, define phase bins,
    amplgmean   = zeros(1,length(phasebins)-1) # preallocate 
    amphgmean   = zeros(1,length(phasebins)-1) # preallocate 
    phithmean   = zeros(1,length(phasebins)-1) # preallocate

    for iBin = 1:length(phasebins)-1
        phLow = phasebins[iBin]   # lower limit
        phHi  = phasebins[iBin+1] # upper limit
        ind   = findall(phLow .<= phith .< phHi)   # find phases falling in bin
        amplgmean[iBin] = mean(amplg[ind])         # compute mean amp
        amphgmean[iBin] = mean(amphg[ind])         # compute mean amp
        phithmean[iBin] = mean([phLow, phHi])      # save center phase 
    end

    hlg = maximum(amplgmean) - minimum(amplgmean) # diff between min and max
    hhg = maximum(amphgmean) - minimum(amphgmean)

    if takepic == 1
        ## plot CFC phase-amplitude
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
        savefig("CFC.png")
    end

    return hlg, hhg 

end