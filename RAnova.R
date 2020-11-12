setwd("D:/DynamicCSDjl")

# load full dataset in
peakdata = read.csv(file = 'AVRECPeakAM.csv')

# cut down to one stim freq on the first stim response
Dat = subset(peakdata, ClickFreq == 2 & OrderofClick ==1)

res.aov3 <- aov(PeakAmp ~ Group * Measurement, data = Dat)
summary(res.aov3)
