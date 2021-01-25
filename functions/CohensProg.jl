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
                # hack the auto alphabatization of the x axis and take only pre and post 1 
                if stimtype[iSti] == "CL"
                    statout.Measurement[statout.Measurement .== "preCL_1"] .= "1PreCL"
                    statout = filter(row -> (row.Measurement == "CL_1" || row.Measurement .== "1PreCL"), statout)
                elseif stimtype[iSti] == "AM"
                    statout.Measurement[statout.Measurement .== "preAM_1"] .= "1PreAM"
                    statout = filter(row -> (row.Measurement == "AM_1" || row.Measurement .== "1PreAM"), statout)
                end

                for iLay = 1:length(layers)
                # pull out current layer
                curLay = statout[statout[:,:Layer] .== layers[iLay],:]

                    if layers[iLay] == "All"
                        Title = "Effect sizes at " * peaks[iPea] * " peak of " * freqtype[iFrq] * " " * stimtype[iSti] * " Avrec"
                    else
                        Title = "Effect sizes at " * peaks[iPea] * " peak of " * freqtype[iFrq] * " " * stimtype[iSti] * " " * layers[iLay]
                    end

                cohensplot = @df curLay groupedbar(:Comparison, :CDRMS, group = :Measurement, bar_position = :dodge, palette = :darktest, ylims=(0,1), title=Title, xlab = "Measurement", ylab = "Cohen's d");

                hline!([0.2,0.5,0.8])

                name = joinpath(figs,foldername,Title) * savetype
                savefig(cohensplot, name)

                #another way to make this figure: 
                cvtpre, cvtpost = curLay[1,9], curLay[2,9]
                cvvpre, cvvpost = curLay[3,9], curLay[4,9]
                tvvpre, tvvpost = curLay[5,9], curLay[6,9]

                cohensplot2 = plot([cvtpre,cvtpost], xlab="Measurement", ylab="Cohen's d", ylims=(0,1), xlims=(0.8,2.2), xticks=1:1:2, markershape=:circle,markersize=10,linewidth=6,label="C vs T",palette=:darktest,size=(200,400))
                plot!([cvvpre, cvvpost],markershape=:circle,markersize=10,linewidth=6,label="C vs V")
                plot!([tvvpre, tvvpost],markershape=:circle,markersize=10,linewidth=6,label="T vs V")
                hline!([0.2,0.5,0.8],label="Effect sizes")

                name = joinpath(figs,foldername,Title) * " line" * savetype
                savefig(cohensplot2, name)
                
                end
            end
        end
    end
end