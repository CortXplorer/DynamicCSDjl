## Set working directory and create necessary file paths
setwd("D:/DynamicCSDjl")
home = getwd()
datapath = file.path(home,"Data")

sink(file=file.path(datapath,"AvrecPeakStats","ANOVASummary.txt"), type=c("output")) #save all output to txt file!
stimtype = list("CL","AM")
comptype = list("KIT vs KIC","KIT vs KIV","KIC vs KIV")

for (iStim in 1:length(stimtype)) {
  print(paste0("=========================== ANOVAs FOR STIM TYPE ", stimtype[iStim], " ==========================="))
  # load full dataset in of appropriate stim type
  peakdata = read.csv(file = file.path(datapath,paste0("AVRECPeak",stimtype[1],".csv")))
  # pull out stim frequencies present to loop through
  stimfreq = unique(peakdata[c("ClickFreq")])
  
  for (iFreq in 1:nrow(stimfreq)) {
    
    print(paste0("===================== STIMULUS FREQUENCY ", stimfreq[1,iFreq], " Hz ====================="))
    
    # cut down to one stim frequency 
    Dat = subset(peakdata, ClickFreq == stimfreq[1,iFreq])
    
    # loop through each response per stimuli 
    for (iOrd in 1:stimfreq[1,iFreq]) {
      
      print(paste0("=============== RESPONSE ", iOrd, " OF ", ordernum, " ==============="))
      Dat = subset(Dat, OrderofClick == iOrd)
      
      for (iCmp in 1:length(comptype)) {
        
        if (comptype[comptype[iCmp]] == "KIT vs KIC") {
          subDat = subset(Dat, !(Group == "KIV"))
        } else if (comptype[comptype[iCmp]] == "KIT vs KIV") {
          subDat = subset(Dat, !(Group == "KIC"))
        } else if (comptype[comptype[iCmp]] == "KIC vs KIV") {
          subDat = subset(Dat, !(Group == "KIT"))
        }
        
        print(paste0("========= PEAK AMPLITUDE ", comptype[iCmp], " =========")) 
        res.aov3 = aov(PeakAmp ~ Group * Measurement, data = subDat)
        summary(res.aov3)
        
        print(paste0("========= PEAK LATITUDE ", comptype[iCmp], " =========")) 
        res.aov3 = aov(PeakLat ~ Group * Measurement, data = subDat)
        summary(res.aov3)
        
        print(paste0("========= ROOT MEAN SQUARE ", comptype[iCmp], " =========")) 
        res.aov3 = aov(RMS ~ Group * Measurement, data = subDat)
        summary(res.aov3)
      
      }
    } # order of click
  } # which stim frequency
} # which stim type, CL or AM

sink() #restore output to console and finish using txt file