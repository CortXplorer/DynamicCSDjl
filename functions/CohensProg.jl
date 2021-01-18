function CohensProg(figs, data, freqtype, stimtype, savetype)
    
    statsf  = joinpath(data,"AvrecPeakStats")
    foldername = "CohensDPlots"
    if !isdir(joinpath(figs,foldername))
        mkdir(joinpath(figs,foldername))
    end

    layers = ["All" "I_II" "IV" "V" "VI"]
    peaks = ["1st" "2nd"]

    for iPea = 1:length(peaks)
        for iFrq = 1:length(freqtype)
            for iSti = 1:length(stimtype)

                # load stats output data
                FileName = peaks[iPea] * "_BetweenGroups_" * freqtype[iFrq] * "_" * stimtype[iSti] * "_ST.csv"
                statout = CSV.File(joinpath(statsf,FileName)) |> DataFrame
                # absolute the CD column for RMS
                statout.CDRMS = abs.(statout.CDRMS)
                # hack the auto alphabatization of the x axis 
                statout.Measurement[statout.Measurement .== "preCL_1"] .= "1PreCL"
                statout.Measurement[statout.Measurement .== "preAM_1"] .= "1PreAM"
                for iLay = 1:length(layers)
                # pull out current layer
                curLay = statout[statout[:,:Layer] .== layers[iLay],:]

                    if layers[iLay] == "All"
                        Title = "Effect sizeds at " * peaks[iPea] * " peak of " * freqtype[iFrq] * " " * stimtype[iSti] * " Avrec"
                    else
                        Title = "Effect sizeds at " * peaks[iPea] * " peak of " * freqtype[iFrq] * " " * stimtype[iSti] * " " * layers[iLay]
                    end

                cohensplot = @df curLay groupedbar(:Measurement, :CDRMS, group = :Comparison, bar_position = :dodge, palette = :darktest, ylims=(0,1), title=Title, xlab = "Measurement", ylab = "Cohen's d");

                hline!([0.2,0.5,0.8])

                name = joinpath(figs,foldername,Title) * savetype
                savefig(cohensplot, name)
                end
            end
        end
    end
end