using PyCall
using CSV, DataFrames

home    = @__DIR__

stimtype = ["CL" "AM"]

sm = pyimport("statsmodels.api") #pyimport_conda("statsmodels.api","statsmodels")
ols = pyimport("statsmodels.formula.api")
pd = pyimport("pandas")

PeakDatTA = pd.read_csv("AVRECPeakCL.csv")

data = PeakDatTA[boop]
data = PeakDatTA[["Group", "Measurement","PeakAmp"]]


# This is 100% not working and I'm only keeping it in case R fails me completely now

# seperate stim frequencies
data = PeakDatTA[PeakDatTA[!,:ClickFreq] .== 2,:]
data = data[data[!,:OrderofClick] .== 1,:]
data = data[:,[:Group,:Measurement,:PeakAmp]]
#data = filter(row -> ! isnan(row.PeakAmp), data)
Dmodel = sm.formula.ols("conformity ~ C(Group,Sum) * C(Measurement,Sum)",data=data).fit()



for iTyp = 1:length(stimtype)
    # Load in data from matlab table csv file which contains 2 and 5 hz peak amp and latency. 
    if stimtype[iTyp] == "CL"
        PeakDatTA = CSV.File("AVRECPeakCL.csv") |> DataFrame # trial average
        #PeakDatST = CSV.File("AVRECPeakCLST.csv") |> DataFrame # single trial
        println("Click Trains")
    elseif stimtype[iTyp] == "AM"
        PeakDatTA = CSV.File("AVRECPeakAM.csv") |> DataFrame # trial average
        #PeakDatST = CSV.File("AVRECPeakAMST.csv") |> DataFrame # single trial
        println("Amplitude Modulation")
    end


end




