

s1BL <- matrix(0, nrow= 250, ncol= 120)
s1BS <- matrix(0, nrow= 250, ncol= 120)
s1TL <- matrix(0, nrow= 250, ncol= 120)
s1TS <- matrix(0, nrow= 250, ncol= 120)
s2BL <- matrix(0, nrow= 250, ncol= 120)
s2BS <- matrix(0, nrow= 250, ncol= 120)
s2TL <- matrix(0, nrow= 250, ncol= 120)
s2TS <- matrix(0, nrow= 250, ncol= 120)
lBS <- matrix(0, nrow= 250, ncol= 120)
lBL <- matrix(0, nrow= 250, ncol= 120)
lTS <- matrix(0, nrow= 250, ncol= 120)
lTL <- matrix(0, nrow= 250, ncol= 120)
infBS <- matrix(0, nrow= 250, ncol= 120)
infBL <- matrix(0, nrow= 250, ncol= 120)
infTS <- matrix(0, nrow= 250, ncol= 120)
infTL <- matrix(0, nrow= 250, ncol= 120)
tBS <- matrix(0, nrow= 250, ncol= 120)
tBL <- matrix(0, nrow= 250, ncol= 120)
tTS <- matrix(0, nrow= 250, ncol= 120)
tTL <- matrix(0, nrow= 250, ncol= 120)
eBS <- matrix(0, nrow= 250, ncol= 120)
eBL <- matrix(0, nrow= 250, ncol= 120)
eTS <- matrix(0, nrow= 250, ncol= 120)
eTL <- matrix(0, nrow= 250, ncol= 120)
p_iBS <- matrix(0, nrow= 250, ncol= 120)
p_iBL <- matrix(0, nrow= 250, ncol= 120)
p_iTS <- matrix(0, nrow= 250, ncol= 120)
p_iTL <- matrix(0, nrow= 250, ncol= 120)
intrcptBS <- matrix(0, nrow= 250, ncol= 120)
intrcptBL <- matrix(0, nrow= 250, ncol= 120)
intrcptTS <- matrix(0, nrow= 250, ncol= 120)
intrcptTL <- matrix(0, nrow= 250, ncol= 120)
runoffBS <- matrix(0, nrow= 250, ncol= 120)
runoffBL <- matrix(0, nrow= 250, ncol= 120)
runoffTS <- matrix(0, nrow= 250, ncol= 120)
runoffTL <- matrix(0, nrow= 250, ncol= 120)
down1BS <- matrix(0, nrow= 250, ncol= 120) 
down2BS <- matrix(0, nrow= 250, ncol= 120)
down1BL <- matrix(0, nrow= 250, ncol= 120) 
down2BL <- matrix(0, nrow= 250, ncol= 120)  
down1TS <- matrix(0, nrow= 250, ncol= 120) 
down2TS <- matrix(0, nrow= 250, ncol= 120)
down1TL <- matrix(0, nrow= 250, ncol= 120) 
down2TL <- matrix(0, nrow= 250, ncol= 120) 
crossing1BS <- matrix(0, nrow= 250, ncol= 119)
crossing2BS <- matrix(0, nrow= 250, ncol= 119)
crossing1BL <- matrix(0, nrow= 250, ncol= 119)
crossing2BL <- matrix(0, nrow= 250, ncol= 119)
crossing1TS <- matrix(0, nrow= 250, ncol= 119)
crossing2TS <- matrix(0, nrow= 250, ncol= 119)
crossing1TL <- matrix(0, nrow= 250, ncol= 119)
crossing2TL <- matrix(0, nrow= 250, ncol= 119)

valueOfInterstormB= matrix(0, nrow= 250, ncol= 120)
valueOfInterstormT= matrix(0, nrow= 250, ncol= 120)

valueOfDurationB= matrix(0, nrow= 250, ncol= 120)
valueOfDurationT= matrix(0, nrow= 250, ncol= 120)

valueOfIntensityB= matrix(0, nrow= 250, ncol= 120)
valueOfIntensityT= matrix(0, nrow= 250, ncol= 120)

inputBS<- rep(0, length(250))
inputBL<- rep(0, length(250))
inputTS<- rep(0, length(250))
inputTL<- rep(0, length(250))

lwbBS<- rep(0, length(250))
lwbBL<- rep(0, length(250))
lwbTS<- rep(0, length(250))
lwbTL<- rep(0, length(250))

twbBS<- rep(0, length(250))
twbBL<- rep(0, length(250))
twbTS<- rep(0, length(250))
twbTL<- rep(0, length(250))

ewbBS<- rep(0, length(250))
ewbBL<- rep(0, length(250))
ewbTS<- rep(0, length(250))
ewbTL<- rep(0, length(250))

iwbBS<- rep(0, length(250))
iwbBL<- rep(0, length(250))
iwbTS<- rep(0, length(250))
iwbTL<- rep(0, length(250))

rwbBS<- rep(0, length(250))
rwbBL<- rep(0, length(250))
rwbTS<- rep(0, length(250))
rwbTL<- rep(0, length(250))

numdown1BS<- rep(0, length(250))
numdown1BL<- rep(0, length(250))
numdown1TS<- rep(0, length(250))
numdown1TL<- rep(0, length(250))


numdown2BS<- rep(0, length(250))
numdown2BL<- rep(0, length(250))
numdown2TS<- rep(0, length(250))
numdown2TL<- rep(0, length(250))

T_bar1BS<- rep(0, length(250))
T_bar1BL<- rep(0, length(250))
T_bar1TS<- rep(0, length(250))
T_bar1TL<- rep(0, length(250))

T_bar2BS<- rep(0, length(250))
T_bar2BL<- rep(0, length(250))
T_bar2TS<- rep(0, length(250))
T_bar2TL<- rep(0, length(250))

#Climate Data
#mtr-duration, mtb-interstorm, mi-intensity
#Boston 
bmtr= 4
bmtb= 68
bmi= 1.691

#Tuscan
tmtr= 44
tmtb= 129
tmi= 0.814

num_days=120

#ETmax
ETmaxB=4
ETmaxT=7
BmaxE= 0.32*ETmaxB
BmaxT= 0.68*ETmaxB
TmaxE= 0.32*ETmaxT
TmaxT= 0.68*ETmaxT

#Sand
nS= 0.42
sfcS= 0.29
saS= 0.105
shS= 0.02
swS= 0.036
bS= 9
ksatS= 1098 #mm/day


#Loam
nL= 0.468
sfcL= 0.592
saL= 0.35
shL= 0.137
swL= 0.165
bL= 9
ksatL= 330 #mm/day

#Vegetation
incpt= 2

zr= 500 
init_s= 0.25


#Rainfall Stimulation
startB <- matrix(, nrow= 250, ncol= 120)
endB <- matrix(, nrow= 250, ncol= 120)
depthB <- matrix(, nrow= 250, ncol= 120)

startT <- matrix(, nrow= 250, ncol= 120)
endT <- matrix(, nrow= 250, ncol= 120)
depthT <- matrix(, nrow= 250, ncol= 120)

PB <- matrix(0, nrow=250, ncol=120)
PT <- matrix(0, nrow=250, ncol=120)#pre-specify daily precipitation variable


for (j in 1:250){

  
  #Random Numbers
  interstorm= runif(num_days)
  valueOfInterstormB[j, ]= -bmtb*log(1-interstorm)
  valueOfInterstormT[j, ]= -tmtb*log(1-interstorm)
  

  duration= runif(num_days)
  valueOfDurationB[j, ]= -bmtr*log(1-duration)
  valueOfDurationT[j, ]= -tmtr*log(1-duration)
  

  intensity= runif(num_days)
  valueOfIntensityB[j, ]= -bmi*log(1-intensity)
  valueOfIntensityT[j, ]= -tmi*log(1-intensity)
  

  
  
  for(i in 1:num_days){
    if(i==1){
      startB[j, i]= valueOfInterstormB[j, i]
      startT[j, i]= valueOfInterstormT[j, i]
    }else{
      startB[j, i]= endB[j, i-1] + valueOfInterstormB[j, i]
      startT[j, i]= endT[j, i-1] + valueOfInterstormT[j, i]
    }
    
    endB[j, i]= startB[j, i] + valueOfDurationB[j, i]
    endT[j, i]= startT[j, i] + valueOfDurationT[j, i]
    
  }
  
  for(i in 1:num_days){
    startB[j, i]=ceiling(startB[j, i]/24)
    endB[j, i]=ceiling(endB[j, i]/24)
    
    startT[j, i]=ceiling(startT[j, i]/24)
    endT[j, i]=ceiling(endT[j, i]/24)
  }
  
  end_dayB=endB[j, ]
  end_dayT=endT[j, ]
  
  depthB= valueOfDurationB[j, ]* valueOfIntensityB[j, ]
  depthT= valueOfDurationT[j, ]* valueOfIntensityT[j, ]
  
  day <- 1:num_days    #days in 120-day time series
  

  
  for (i in 1:length(day)) {
    
    #find any days with rain storms ending on the given day
    II <- which(end_dayB == day[i]) 
    HH <- which(end_dayT == day[i])
    
    if (length(II) >= 1) {
      PB[j, i] <- sum(depthB[II])    #calculate total rain depth on the given day, and set value
    }
    if (length(HH) >= 1) {
      PT[j, i] <- sum(depthT[HH])    #calculate total rain depth on the given day, and set value
    }
  }
  #plot(day, PB[j], type= "h", xlab="Day", ylab= "Precipitation (mm)", main= "Rainfall in Boston")
  #plot(day, PT[j], type= "h", xlab="Day", ylab= "Precipitation (mm)", main= "Rainfall in Tuscon")
  
  
  for (i in 1:num_days){
    #Precipitaion-Interception
    p_iBS[j, i]=PB[j, i]- incpt
    p_iBL[j, i]=PB[j, i]- incpt
    p_iTS[j, i]=PT[j, i]- incpt
    p_iTL[j, i]=PT[j, i]- incpt
    if(p_iBS[j, i]<0){
      p_iBS[j, i]=0
    }
    if(p_iBL[j, i]<0){
      p_iBL[j, i]=0
    }
    if(p_iTS[j, i]<0){
      p_iTS[j, i]=0
    }
    if(p_iTL[j, i]<0){
      p_iTL[j, i]=0
    }
    
    #Infiltration
    if(i==1){
      infBS[j, i]= min(p_iBS[j, i], nS*zr*(1-init_s))
      infBL[j, i]= min(p_iBL[j, i], nL*zr*(1-init_s))
      infTS[j, i]= min(p_iTS[j, i], nS*zr*(1-init_s))
      infTL[j, i]= min(p_iTL[j, i], nL*zr*(1-init_s))
    }else{
      infBS[j, i]= min(p_iBS[j, i], nS*zr*(1-s2BS[j, i-1]))
      infBL[j, i]= min(p_iBL[j, i], nL*zr*(1-s2BL[j, i-1]))
      infTS[j, i]= min(p_iTS[j, i], nS*zr*(1-s2TS[j, i-1]))
      infTL[j, i]= min(p_iTL[j, i], nL*zr*(1-s2TL[j, i-1]))
    }
    
    #Saturation1
    if(i==1){
      s1BS[j, i]= (infBS[j, i]+(nS*zr*init_s))/(nS*zr)
      s1BL[j, i]= (infBL[j, i]+(nL*zr*init_s))/(nL*zr)
      s1TS[j, i]= (infTS[j, i]+(nS*zr*init_s))/(nS*zr)
      s1TL[j, i]= (infTL[j, i]+(nL*zr*init_s))/(nL*zr)
    }else{
      s1BS[j, i]= (infBS[j, i]+(nS*zr*s2BS[j, i-1]))/(nS*zr)
      s1BL[j, i]= (infBL[j, i]+(nL*zr*s2BL[j, i-1]))/(nL*zr)
      s1TS[j, i]= (infTS[j, i]+(nS*zr*s2TS[j, i-1]))/(nS*zr)
      s1TL[j, i]= (infTL[j, i]+(nL*zr*s2TL[j, i-1]))/(nL*zr)
    }
    #Leakage
    if(s1BS[j, i]<sfcS){
      lBS[j, i]=0
    }else{
      lBS[j, i]= min(ksatS*((exp(bS*(s1BS[j, i]-sfcS))-1)/(exp(bS*(1-sfcS))-1)), nS*zr*(s1BS[j, i]-sfcS))
    }
    if(s1BL[j, i]<sfcL){
      lBL[j, i]=0
    }else{
      lBL[j, i]= min(ksatL*((exp(bL*(s1BL[j, i]-sfcL))-1)/(exp(bL*(1-sfcL))-1)), nL*zr*(s1BL[j, i]-sfcL))
    }
    if(s1TS[j, i]<sfcS){
      lTS[j, i]=0
    }else{
      lTS[j, i]= min(ksatS*((exp(bS*(s1TS[j, i]-sfcS))-1)/(exp(bS*(1-sfcS))-1)), nS*zr*(s1TS[j, i]-sfcS))
    }
    if(s1TL[j, i]<sfcL){
      lTL[j, i]=0
    }else{
      lTL[j, i]= min(ksatL*((exp(bL*(s1TL[j, i]-sfcL))-1)/(exp(bL*(1-sfcL))-1)), nL*zr*(s1TL[j, i]-sfcL))
    }
    
    if(s1BS[j, i]<=swS){
      tBS[j, i]=0
    }
    if(s1BL[j, i]<=swL){
      tBL[j, i]=0
    }
    if(s1TS[j, i]<=swS){
      tTS[j, i]=0
    }
    if(s1TL[j, i]<=swL){
      tTL[j, i]=0
    }
    
    if(swS<s1BS[j, i] & s1BS[j, i]<saS){
      tBS[j, i]= ((s1BS[j, i]-swS)/(saS-swS))*BmaxT
    }
    if(swL<s1BL[j, i] & s1BL[j, i]<saL){
      tBL[j, i]= ((s1BL[j, i]-swL)/(saL-swL))*BmaxT
    }
    if(swS<s1TS[j, i] & s1TS[j, i]<saS){
      tTS[j, i]= ((s1TS[j, i]-swS)/(saS-swS))*TmaxT
    }
    if(swL<s1TL[j, i] & s1TL[j, i]<saL){
      tTL[j, i]= ((s1TL[j, i]-swL)/(saL-swL))*TmaxT
    }
    
    if(s1BS[j, i]>=saS){
      tBS[j, i]=BmaxT
    }
    if(s1BL[j, i]>=saL){
      tBL[j, i]=BmaxT
    }
    if(s1TS[j, i]>=saS){
      tTS[j, i]=TmaxT
    }
    if(s1TL[j, i]>=saL){
      tTL[j, i]=TmaxT
    }
    
    if(s1BS[j, i]<=shS){
      eBS[j, i]=0
    }
    if(s1TS[j, i]<=shS){
      eTS[j, i]=0
    }
    if(s1BL[j, i]<=shL){
      eBL[j, i]=0
    }
    if(s1TL[j, i]<=shL){
      eTL[j, i]=0
    }
    if(shS<s1BS[j, i] & s1BS[j, i]<saS){
      eBS[j, i]= ((s1BS[j, i]-shS)/(saS-shS))*BmaxE
    }
    if(shL<s1BL[j, i] & s1BL[j, i]<saL){
      eBL[j, i]= ((s1BL[j, i]-shL)/(saL-shL))*BmaxE
    }
    if(shS<s1TS[j, i] & s1TS[j, i]<saS){
      eTS[j, i]= ((s1TS[j, i]-shS)/(saS-shS))*TmaxE
    }
    if(shL<s1TL[j, i] & s1TL[j, i]<saL){
      eTL[j, i]= ((s1TL[j, i]-shL)/(saL-shL))*TmaxE
    }
    if(s1BS[j, i]>=saS){
      eBS[j, i]= BmaxE
    }
    if(s1BL[j, i]>=saL){
      eBL[j, i]= BmaxE
    }
    if(s1TS[j, i]>=saS){
      eTS[j, i]= TmaxE
    }
    if(s1TL[j, i]>=saL){
      eTL[j, i]= TmaxE
    }
    
    if(PB[j, i]>0){
      #Transpiration
      tBS[j, i]=0
      tBL[j, i]=0
      #Evaporation
      eBS[j, i]=0
      eBL[j, i]=0
    }
    if(PT[j, i]>0){
      #Transpiration
      tTS[j, i]=0
      tTL[j, i]=0
      #Evaporation
      eTS[j, i]=0
      eTL[j, i]=0
    }
    
    #Saturation2
    s2BS[j, i]= ((nS*zr*s1BS[j, i])-lBS[j, i]-tBS[j, i]-eBS[j, i])/(nS*zr)
    s2BL[j, i]= ((nL*zr*s1BL[j, i])-lBL[j, i]-tBL[j, i]-eBL[j, i])/(nL*zr)
    s2TS[j, i]= ((nS*zr*s1TS[j, i])-lTS[j, i]-tTS[j, i]-eTS[j, i])/(nS*zr)
    s2TL[j, i]= ((nL*zr*s1TL[j, i])-lTL[j, i]-tTL[j, i]-eTL[j, i])/(nL*zr)
  }
  
  
#   plot(day, s2BS, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Boston and Sandy")
#   abline(h=swS, col="red", lty=2)
#   abline(h=saS, col="blue", lty=2)
#   abline(h=sfcS, col="green", lty=2)
#   legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")
#   
#   plot(day, s2BL, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Boston and Silty Loam")
#   abline(h=swL, col="red", lty=2)
#   abline(h=saL, col="blue", lty=2)
#   abline(h=sfcL, col="green", lty=2)
#   legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")
#   
#   plot(day, s2TS, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Tucson and Sandy")
#   abline(h=swS, col="red", lty=2)
#   abline(h=saS, col="blue", lty=2)
#   abline(h=sfcS, col="green", lty=2)
#   legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")
#   
#   plot(day, s2TL, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Tucson and Silty Loam")
#   abline(h=swL, col="red", lty=2)
#   abline(h=saL, col="blue", lty=2)
#   abline(h=sfcL, col="green", lty=2)
#   legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")
#   
#   
  
  #Interception
  intrcptBS[j, ]= PB[j, ]-p_iBS[j, ]
  intrcptBL[j, ]= PB[j, ]-p_iBL[j, ]
  intrcptTS[j, ]= PT[j, ]-p_iTS[j, ]
  intrcptTL[j, ]= PT[j, ]-p_iTL[j, ]
  
  #Runoff
  runoffBS[j, ]= p_iBS[j, ]-infBS[j, ]
  runoffBL[j, ]= p_iBL[j, ]-infBL[j, ]
  runoffTS[j, ]= p_iTS[j, ]-infTS[j, ]
  runoffTL[j, ]= p_iTL[j, ]-infTL[j, ]
  
  #input
  inputBS[j]= sum(PB[j, ])-nS*zr*(s2BS[j, num_days]-init_s)
  inputBL[j]= sum(PB[j, ])-nL*zr*(s2BL[j, num_days]-init_s)
  inputTS[j]= sum(PT[j, ])-nS*zr*(s2TS[j, num_days]-init_s)
  inputTL[j]= sum(PT[j, ])-nL*zr*(s2TL[j, num_days]-init_s)
  
  #water balance: leakage
  lwbBS[j]= sum(lBS[j, ])/inputBS[j]
  lwbBL[j]= sum(lBL[j, ])/inputBL[j]
  lwbTS[j]= sum(lTS[j, ])/inputTS[j]
  lwbTL[j]= sum(lTL[j, ])/inputTL[j]
  
  #water balance: transpirtion
  twbBS[j]= sum(tBS[j, ])/inputBS[j]
  twbBL[j]= sum(tBL[j, ])/inputBL[j]
  twbTS[j]= sum(tTS[j, ])/inputTS[j]
  twbTL[j]= sum(tTL[j, ])/inputTL[j]
  
  #water balance: evaporation
  ewbBS[j]= sum(eBS[j, ])/inputBS[j]
  ewbBL[j]= sum(eBL[j, ])/inputBL[j]
  ewbTS[j]= sum(eTS[j, ])/inputTS[j]
  ewbTL[j]= sum(eTL[j, ])/inputTL[j]
  
  #water balance: interception
  iwbBS[j]= sum(intrcptBS[j, ])/inputBS[j]
  iwbBL[j]= sum(intrcptBL[j, ])/inputBL[j]
  iwbTS[j]= sum(intrcptTS[j, ])/inputTS[j]
  iwbTL[j]= sum(intrcptTL[j, ])/inputTL[j]
  
  #water balance: runoff
  rwbBS[j]= sum(runoffBS[j, ])/inputBS[j]
  rwbBL[j]= sum(runoffBL[j, ])/inputBL[j]
  rwbTS[j]= sum(runoffTS[j, ])/inputTS[j]
  rwbTL[j]= sum(runoffTL[j, ])/inputTL[j]
  
  lwbBS[j]+twbBS[j]+ewbBS[j]+iwbBS[j]+rwbBS[j]
  lwbBL[j]+twbBL[j]+ewbBL[j]+iwbBL[j]+rwbBL[j]
  lwbTS[j]+twbTS[j]+ewbTS[j]+iwbTS[j]+rwbTS[j]
  lwbTL[j]+twbTL[j]+ewbTL[j]+iwbTL[j]+rwbTL[j]
  
  
#   counts <- matrix(c(lwbBS, lwbBL, lwbTS, lwbTL, twbBS, twbBL, twbTS, twbTL, ewbBS, ewbBL, ewbTS, ewbTL, iwbBS, iwbBL, iwbTS, iwbTL, rwbBS, rwbBL, rwbTS, rwbTL), nrow=5, byrow=TRUE)      
#   colnames(counts) <- c("BS","BL","TS", "TL")
#   rownames(counts) <- c("l", "t", "e", "i", "r")
#   counts <- as.table(counts)
#   barplot(counts,legend=T,beside=F,main='Water Output Balance', col=c("blue", "red", "green", "coral", "darkcyan"))
  
  
  
  
#   bar1 <- matrix(c(num_down1BS[j], num_down1BL[j], num_down1TS[j], num_down1TL[j], num_down2BS[j], num_down2BL[j], num_down2TS[j], num_down2TL[j]), nrow=2, byrow=TRUE)      
#   colnames(bar1) <- c("BS","BL","TS", "TL")
#   rownames(bar1) <- c("num_down1", "num_down2")
#   bar <- as.table(bar1)
  #barplot(bar1,legend=T,beside=T,main='Number of Downcrossings', col=c("red", "coral"))
  

  
#   bar2 <- matrix(c( T_bar1BS[j], T_bar1BL[j], T_bar1TS[j], T_bar1TL[j], T_bar2BS[j], T_bar2BL[j], T_bar2TS[j], T_bar2TL[j]), nrow=2, byrow=TRUE)      
#   colnames(bar2) <- c("BS","BL","TS", "TL")
#   rownames(bar2) <- c("T_bar1", "T_bar2")
#   bar <- as.table(bar2)
  #barplot(bar2,legend=T,beside=T,main='Duration of Excursion', col=c("blue", "cyan"))
  
  
  
  
  
  
}

(sum(lwbBS)+sum(twbBS)+sum(ewbBS)+sum(iwbBS)+sum(rwbBS))/250
(sum(lwbBL)+sum(twbBL)+sum(ewbBL)+sum(iwbBL)+sum(rwbBL))/250
(sum(lwbTS)+sum(twbTS)+sum(ewbTS)+sum(iwbTS)+sum(rwbTS))/250
(sum(lwbTL)+sum(twbTL)+sum(ewbTL)+sum(iwbTL)+sum(rwbTL))/250

l_BS= sum(lwbBS)/250
l_BL= sum(lwbBL)/250
l_TS= sum(lwbTS)/250
l_TL= sum(lwbTL)/250

t_BS= sum(twbBS)/250
t_BL= sum(twbBL)/250
t_TS= sum(twbTS)/250
t_TL= sum(twbTL)/250

e_BS= sum(ewbBS)/250
e_BL= sum(ewbBL)/250
e_TS= sum(ewbTS)/250
e_TL= sum(ewbTL)/250

i_BS= sum(iwbBS)/250
i_BL= sum(iwbBL)/250
i_TS= sum(iwbTS)/250
i_TL= sum(iwbTL)/250

r_BS= sum(rwbBS)/250
r_BL= sum(rwbBL)/250
r_TS= sum(rwbTS)/250
r_TL= sum(rwbTL)/250

  counts <- matrix(c(l_BS, l_BL, l_TS, l_TL, t_BS, t_BL, t_TS, t_TL, e_BS, e_BL, e_TS, e_TL, i_BS, i_BL, i_TS, i_TL, r_BS, r_BL, r_TS, r_TL), nrow=5, byrow=TRUE)      
  colnames(counts) <- c("BS","BL","TS", "TL")
  rownames(counts) <- c("l", "t", "e", "i", "r")
  counts <- as.table(counts)
  barplot(counts,legend=T,beside=F,main='Average Water Output Balance', col=c("blue", "red", "green", "coral", "darkcyan"))

