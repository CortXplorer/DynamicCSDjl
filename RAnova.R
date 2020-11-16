## Set working directory and create necessary file paths
setwd("D:/DynamicCSDjl")
home = getwd()
datapath = file.path(home,"Data")

# load full dataset in
peakdata = read.csv(file = file.path(datapath,"AVRECPeakAM.csv"))

# cut down to one stim freq on the first stim response
Dat = subset(peakdata, ClickFreq == 2 & OrderofClick ==1)

res.aov3 <- aov(PeakAmp ~ Group * Measurement, data = Dat)
DatSum = summary(res.aov3)

#method one to save out one at a time
capture.output(DatSum,file="Datsum.txt",append = FALSE,type = c("output","message"),split=FALSE)

# method 2 of creating print lines as titles and data output for all dat to one file
sink("DataSummary.txt", type=c("output"))
print("Hello there boss")
DatSum
sink()
