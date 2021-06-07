
function getAMwave(clickfreq=2,sr=1000)

    t  = 0:1/sr:1 # time course
    wave = sin.((2*π*clickfreq*t).+3π/2)
    return wave

end

function VSstats(data, sr, type, clickfreq)
    # open the dataframe in the workspace
    nametype = "AVRECPeak" * type * "ST.csv"
    PeakDatST = CSV.File(joinpath(data,nametype)) |> DataFrame
    PeakDatST = PeakDatST[PeakDatST[:,:ClickFreq] .== clickfreq,:] 

    # set up DataFrame:
    df = DataFrame(Group = String[], Animal = String[], Layer = String[], Measurement = String[], VectorStrenght = Float64[], P = Float64[], Z = Float64[], MeanPhase = Float64[])

    # loop through parameters we're including
    for iGr = unique(PeakDatST[:,:Group])
        GrDat = PeakDatST[PeakDatST[:,:Group] .== iGr,:]
        for iAn = unique(GrDat[:,:Animal])
            AnDat = GrDat[GrDat[:,:Animal] .== iAn,:]
            for iLay = unique(AnDat[:,:Layer])
                LaDat = AnDat[AnDat[:,:Layer] .== iLay,:]
                for iMea = unique(LaDat[:,:Measurement])
                    CurDat = LaDat[LaDat[:,:Measurement] .== iMea,:]

                    # our final list of latencies:
                    CurDat      = filter(row -> ! isnan(row.PeakLat), CurDat)
                    peaklat     = Int.(CurDat[:,:PeakLat])
                    orderlist   = Int.(CurDat[:,:OrderofClick])
                    # get the right timing based on stim presentation and remove responses before 100 ms to avoid onset synchronization
                    peaklat_cor = ((orderlist .- 1) .* Int(1000/clickfreq)) .+ peaklat
                    peaklat_cor = peaklat_cor[peaklat_cor .> 100]

                    # generate the AM wave used and get instantaneous phase
                    wave = getAMwave(clickfreq,sr) # currently below
                    wavephase = angle.(hilbert(wave)) # Compute phase of wave

                    # get the phase at each time point of signal peak amplitude
                    orientation = wavephase[peaklat_cor]

                    ### Rayleigh Test of circular uniformity ### 
                    # ______________________________________________________________________
                    # break vector into components i and j to compute sums, averages, and then the resultant vector 
                    i = cos.(orientation) # get component i 
                    j = sin.(orientation)
                    sumi = sum(i) # 
                    sumj = sum(j)
                    avgi = sumi/length(i) # average of component i
                    avgj = sumj/length(j) # average of component j

                    vstrength = sqrt((avgi^2)+(avgj^2)) # find the vector based on i and j components (pythagoreon)

                    # zscore and subsequent p value:
                    Z = vstrength^2 * length(orientation)
                    p = exp(-Z) 

                    # we can plot that as a simple histogram: 
                    # histogram(orientation,xlims=(-π,π),xticks=-π:4/π:π,label="moderator phase at peak latency")
                    # ______________________________________________________________________

                    meanphase = atan(abs(avgj)/abs(avgi)) # find the mean phase in radian; 
                    # this gives the corresponding angle upper right quadrant. we need to translate it to the correct quadrant based on the polar graph and signs of i and j
                    if avgi > 0 && avgj > 0     # upper right target
                        meanphase = 0 + meanphase
                    elseif avgi > 0 && avgj < 0 # lower right target
                        meanphase = 2π - meanphase
                    elseif avgi < 0 && avgj < 0 # lower left target
                        meanphase = π + meanphase
                    elseif avgi < 0 && avgj > 0 # upper left target
                        meanphase = π - meanphase
                    end

                    push!(df, [iGr  iAn iLay iMea vstrength p Z meanphase])
                end
            end
        end
    end
    foldername = "VSoutput"
    if !isdir(joinpath(data,foldername))
        mkdir(joinpath(data,foldername))
    end
    title = type * "_" * string(clickfreq) * "_VectorStrength_Stats" 
    name = joinpath(data,foldername,title) * ".csv"
    CSV.write(name, df)
end