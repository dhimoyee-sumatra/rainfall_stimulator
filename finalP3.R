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

#Random Numbers
set.seed(48252)
interstorm= runif(num_days)
valueOfInterstormB= -bmtb*log(1-interstorm)
valueOfInterstormT= -tmtb*log(1-interstorm)

set.seed(39616)
duration= runif(num_days)
valueOfDurationB= -bmtr*log(1-duration)
valueOfDurationT= -tmtr*log(1-duration)

set.seed(41106)
intensity= runif(num_days)
valueOfIntensityB= -bmi*log(1-intensity)
valueOfIntensityT= -tmi*log(1-intensity)

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
startB <- rep(NA, length(num_days))
endB <- rep(NA, length(num_days))
depthB <- rep(NA, length(num_days))

startT <- rep(NA, length(num_days))
endT <- rep(NA, length(num_days))
depthT <- rep(NA, length(num_days))


for(i in 1:num_days){
  if(i==1){
    startB[i]= valueOfInterstormB[i]
    startT[i]= valueOfInterstormT[i]
  }else{
    startB[i]= endB[i-1] + valueOfInterstormB[i]
    startT[i]= endT[i-1] + valueOfInterstormT[i]
  }
  
  endB[i]= startB[i] + valueOfDurationB[i]
  endT[i]= startT[i] + valueOfDurationT[i]
  
}

for(i in 1:num_days){
  startB[i]=ceiling(startB[i]/24)
  endB[i]=ceiling(endB[i]/24)
  
  startT[i]=ceiling(startT[i]/24)
  endT[i]=ceiling(endT[i]/24)
}

end_dayB=endB
end_dayT=endT

depthB= valueOfDurationB* valueOfIntensityB
depthT= valueOfDurationT* valueOfIntensityT

day <- 1:num_days    #days in 120-day time series

PB <- rep(0,length(day))
PT <- rep(0,length(day))#pre-specify daily precipitation variable

for (i in 1:length(day)) {
  
  #find any days with rain storms ending on the given day
  II <- which(end_dayB == day[i]) 
  HH <- which(end_dayT == day[i])
  
  if (length(II) >= 1) {
    PB[i] <- sum(depthB[II])    #calculate total rain depth on the given day, and set value
  }
  if (length(HH) >= 1) {
    PT[i] <- sum(depthT[HH])    #calculate total rain depth on the given day, and set value
  }
}
plot(day, PB, type= "h", xlab="Day", ylab= "Precipitation (mm)", main= "Rainfall in Boston")
plot(day, PT, type= "h", xlab="Day", ylab= "Precipitation (mm)", main= "Rainfall in Tuscon")

s1BL <- rep(NA, length(num_days))
s1BS <- rep(NA, length(num_days))
s1TL <- rep(NA, length(num_days))
s1TS <- rep(NA, length(num_days))
s2BL <- rep(NA, length(num_days))
s2BS <- rep(NA, length(num_days))
s2TL <- rep(NA, length(num_days))
s2TS <- rep(NA, length(num_days))
lBS <- rep(NA, length(num_days))
lBL <- rep(NA, length(num_days))
lTS <- rep(NA, length(num_days))
lTL <- rep(NA, length(num_days))
infBS <- rep(NA, length(num_days))
infBL <- rep(NA, length(num_days))
infTS <- rep(NA, length(num_days))
infTL <- rep(NA, length(num_days))
tBS <- rep(NA, length(num_days))
tBL <- rep(NA, length(num_days))
tTS <- rep(NA, length(num_days))
tTL <- rep(NA, length(num_days))
eBS <- rep(NA, length(num_days))
eBL <- rep(NA, length(num_days))
eTS <- rep(NA, length(num_days))
eTL <- rep(NA, length(num_days))
p_iBS <- rep(NA, length(num_days))
p_iBL <- rep(NA, length(num_days))
p_iTS <- rep(NA, length(num_days))
p_iTL <- rep(NA, length(num_days))

for (i in 1:num_days){
  #Precipitaion-Interception
  p_iBS[i]=PB[i]- incpt
  p_iBL[i]=PB[i]- incpt
  p_iTS[i]=PT[i]- incpt
  p_iTL[i]=PT[i]- incpt
  if(p_iBS[i]<0){
    p_iBS[i]=0
  }
  if(p_iBL[i]<0){
    p_iBL[i]=0
  }
  if(p_iTS[i]<0){
    p_iTS[i]=0
  }
  if(p_iTL[i]<0){
    p_iTL[i]=0
  }
  
  #Infiltration
  if(i==1){
    infBS[i]= min(p_iBS[i], nS*zr*(1-init_s))
    infBL[i]= min(p_iBL[i], nL*zr*(1-init_s))
    infTS[i]= min(p_iTS[i], nS*zr*(1-init_s))
    infTL[i]= min(p_iTL[i], nL*zr*(1-init_s))
  }else{
    infBS[i]= min(p_iBS[i], nS*zr*(1-s2BS[i-1]))
    infBL[i]= min(p_iBL[i], nL*zr*(1-s2BL[i-1]))
    infTS[i]= min(p_iTS[i], nS*zr*(1-s2TS[i-1]))
    infTL[i]= min(p_iTL[i], nL*zr*(1-s2TL[i-1]))
  }
  
  #Saturation1
  if(i==1){
    s1BS[i]= (infBS[i]+(nS*zr*init_s))/(nS*zr)
    s1BL[i]= (infBL[i]+(nL*zr*init_s))/(nL*zr)
    s1TS[i]= (infTS[i]+(nS*zr*init_s))/(nS*zr)
    s1TL[i]= (infTL[i]+(nL*zr*init_s))/(nL*zr)
  }else{
    s1BS[i]= (infBS[i]+(nS*zr*s2BS[i-1]))/(nS*zr)
    s1BL[i]= (infBL[i]+(nL*zr*s2BL[i-1]))/(nL*zr)
    s1TS[i]= (infTS[i]+(nS*zr*s2TS[i-1]))/(nS*zr)
    s1TL[i]= (infTL[i]+(nL*zr*s2TL[i-1]))/(nL*zr)
  }
  #Leakage
  if(s1BS[i]<sfcS){
    lBS[i]=0
  }else{
    lBS[i]= min(ksatS*((exp(bS*(s1BS[i]-sfcS))-1)/(exp(bS*(1-sfcS))-1)), nS*zr*(s1BS[i]-sfcS))
  }
  if(s1BL[i]<sfcL){
    lBL[i]=0
  }else{
    lBL[i]= min(ksatL*((exp(bL*(s1BL[i]-sfcL))-1)/(exp(bL*(1-sfcL))-1)), nL*zr*(s1BL[i]-sfcL))
  }
  if(s1TS[i]<sfcS){
    lTS[i]=0
  }else{
    lTS[i]= min(ksatS*((exp(bS*(s1TS[i]-sfcS))-1)/(exp(bS*(1-sfcS))-1)), nS*zr*(s1TS[i]-sfcS))
  }
  if(s1TL[i]<sfcL){
    lTL[i]=0
  }else{
    lTL[i]= min(ksatL*((exp(bL*(s1TL[i]-sfcL))-1)/(exp(bL*(1-sfcL))-1)), nL*zr*(s1TL[i]-sfcL))
  }
  
  if(s1BS[i]<=swS){
    tBS[i]=0
  }
  if(s1BL[i]<=swL){
    tBL[i]=0
  }
  if(s1TS[i]<=swS){
    tTS[i]=0
  }
  if(s1TL[i]<=swL){
    tTL[i]=0
  }
  
  if(swS<s1BS[i] & s1BS[i]<saS){
    tBS[i]= ((s1BS[i]-swS)/(saS-swS))*BmaxT
  }
  if(swL<s1BL[i] & s1BL[i]<saL){
    tBL[i]= ((s1BL[i]-swL)/(saL-swL))*BmaxT
  }
  if(swS<s1TS[i] & s1TS[i]<saS){
    tTS[i]= ((s1TS[i]-swS)/(saS-swS))*TmaxT
  }
  if(swL<s1TL[i] & s1TL[i]<saL){
    tTL[i]= ((s1TL[i]-swL)/(saL-swL))*TmaxT
  }
  
  if(s1BS[i]>=saS){
    tBS[i]=BmaxT
  }
  if(s1BL[i]>=saL){
    tBL[i]=BmaxT
  }
  if(s1TS[i]>=saS){
    tTS[i]=TmaxT
  }
  if(s1TL[i]>=saL){
    tTL[i]=TmaxT
  }
  
  if(s1BS[i]<=shS){
    eBS[i]=0
  }
  if(s1TS[i]<=shS){
    eTS[i]=0
  }
  if(s1BL[i]<=shL){
    eBL[i]=0
  }
  if(s1TL[i]<=shL){
    eTL[i]=0
  }
  if(shS<s1BS[i] & s1BS[i]<saS){
    eBS[i]= ((s1BS[i]-shS)/(saS-shS))*BmaxE
  }
  if(shL<s1BL[i] & s1BL[i]<saL){
    eBL[i]= ((s1BL[i]-shL)/(saL-shL))*BmaxE
  }
  if(shS<s1TS[i] & s1TS[i]<saS){
    eTS[i]= ((s1TS[i]-shS)/(saS-shS))*TmaxE
  }
  if(shL<s1TL[i] & s1TL[i]<saL){
    eTL[i]= ((s1TL[i]-shL)/(saL-shL))*TmaxE
  }
  if(s1BS[i]>=saS){
    eBS[i]= BmaxE
  }
  if(s1BL[i]>=saL){
    eBL[i]= BmaxE
  }
  if(s1TS[i]>=saS){
    eTS[i]= TmaxE
  }
  if(s1TL[i]>=saL){
    eTL[i]= TmaxE
  }
  
  if(PB[i]>0){
    #Transpiration
    tBS[i]=0
    tBL[i]=0
    #Evaporation
    eBS[i]=0
    eBL[i]=0
  }
  if(PT[i]>0){
    #Transpiration
    tTS[i]=0
    tTL[i]=0
    #Evaporation
    eTS[i]=0
    eTL[i]=0
  }
  
  #Saturation2
  s2BS[i]= ((nS*zr*s1BS[i])-lBS[i]-tBS[i]-eBS[i])/(nS*zr)
  s2BL[i]= ((nL*zr*s1BL[i])-lBL[i]-tBL[i]-eBL[i])/(nL*zr)
  s2TS[i]= ((nS*zr*s1TS[i])-lTS[i]-tTS[i]-eTS[i])/(nS*zr)
  s2TL[i]= ((nL*zr*s1TL[i])-lTL[i]-tTL[i]-eTL[i])/(nL*zr)
}


plot(day, s2BS, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Boston and Sandy")
abline(h=swS, col="red", lty=2)
abline(h=saS, col="blue", lty=2)
abline(h=sfcS, col="green", lty=2)
legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")

plot(day, s2BL, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Boston and Silty Loam")
abline(h=swL, col="red", lty=2)
abline(h=saL, col="blue", lty=2)
abline(h=sfcL, col="green", lty=2)
legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")

plot(day, s2TS, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Tucson and Sandy")
abline(h=swS, col="red", lty=2)
abline(h=saS, col="blue", lty=2)
abline(h=sfcS, col="green", lty=2)
legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")

plot(day, s2TL, type= "l", xlab="Day", ylab= "Soil Moisture", main= "Tucson and Silty Loam")
abline(h=swL, col="red", lty=2)
abline(h=saL, col="blue", lty=2)
abline(h=sfcL, col="green", lty=2)
legend(95, 0.33, legend=c("S*", "Sw", "Sfc"),col=c("blue", "red", "green"), lty=1:2, cex=0.6,bty = "n")



#Interception
intrcptBS= PB-p_iBS
intrcptBL= PB-p_iBL
intrcptTS= PT-p_iTS
intrcptTL= PT-p_iTL

#Runoff
runoffBS= p_iBS-infBS
runoffBL= p_iBL-infBL
runoffTS= p_iTS-infTS
runoffTL= p_iTL-infTL

#input
inputBS= sum(PB)-nS*zr*(s2BS[num_days]-init_s)
inputBL= sum(PB)-nL*zr*(s2BL[num_days]-init_s)
inputTS= sum(PT)-nS*zr*(s2TS[num_days]-init_s)
inputTL= sum(PT)-nL*zr*(s2TL[num_days]-init_s)

#water balance: leakage
lwbBS= sum(lBS)/inputBS
lwbBL= sum(lBL)/inputBL
lwbTS= sum(lTS)/inputTS
lwbTL= sum(lTL)/inputTL

#water balance: transpirtion
twbBS= sum(tBS)/inputBS
twbBL= sum(tBL)/inputBL
twbTS= sum(tTS)/inputTS
twbTL= sum(tTL)/inputTL

#water balance: evaporation
ewbBS= sum(eBS)/inputBS
ewbBL= sum(eBL)/inputBL
ewbTS= sum(eTS)/inputTS
ewbTL= sum(eTL)/inputTL

#water balance: interception
iwbBS= sum(intrcptBS)/inputBS
iwbBL= sum(intrcptBL)/inputBL
iwbTS= sum(intrcptTS)/inputTS
iwbTL= sum(intrcptTL)/inputTL

#water balance: runoff
rwbBS= sum(runoffBS)/inputBS
rwbBL= sum(runoffBL)/inputBL
rwbTS= sum(runoffTS)/inputTS
rwbTL= sum(runoffTL)/inputTL

lwbBS+twbBS+ewbBS+iwbBS+rwbBS
lwbBL+twbBL+ewbBL+iwbBL+rwbBL
lwbTS+twbTS+ewbTS+iwbTS+rwbTS
lwbTL+twbTL+ewbTL+iwbTL+rwbTL


counts <- matrix(c(lwbBS, lwbBL, lwbTS, lwbTL, twbBS, twbBL, twbTS, twbTL, ewbBS, ewbBL, ewbTS, ewbTL, iwbBS, iwbBL, iwbTS, iwbTL, rwbBS, rwbBL, rwbTS, rwbTL), nrow=5, byrow=TRUE)      
colnames(counts) <- c("BS","BL","TS", "TL")
rownames(counts) <- c("l", "t", "e", "i", "r")
counts <- as.table(counts)
barplot(counts,legend=T,beside=F,main='Water Output Balance', col=c("blue", "red", "green", "coral", "darkcyan"))



#number of downcrossings
IIBS <- which(s2BS < swS)    #find which timesteps have an S value below S_crit
HHBS <- which(s2BS < saS)

#make a variable that has a value of 1 any time S is below S_crit, and is zero otherwise
down1BS <- rep(0,length(s2BS)) 
down2BS <- rep(0,length(s2BS)) 
down1BS[IIBS] <- 1
down2BS[HHBS] <- 1


#make a variable that has a value of 1 any time S crosses below S_crit, and is zero otherwise
crossing1BS <- rep(0,length(s2BS)-1)
crossing2BS <- rep(0,length(s2BS)-1)
for (i in 2:length(s2BS)) {
  #if S was above S_crit in the previous timestep, but is below S_crit in the current
  #timestep, a downcrossing occurred - set "crossing" = 1 for this timestep
  if (down1BS[i-1] ==0 & down1BS[i] ==1) {  
    crossing1BS[i] <- 1    #downcrossing
  }
  if (down2BS[i-1] ==0 & down2BS[i] ==1) {  
    crossing2BS[i] <- 1    #downcrossing
  }
}

num_down1BS <- sum(crossing1BS)    #total number of downcrossings
num_down2BS <- sum(crossing2BS) 

#mean length of downcrossings (units: days)
T_bar1BS <- sum(down1BS) / num_down1BS
T_bar2BS <- sum(down2BS) / num_down2BS
#number of days that S was below S_crit divided by number of downcrossings



#number of downcrossings
IIBL <- which(s2BL < swL)    #find which timesteps have an S value below S_crit
HHBL <- which(s2BL < saL)

#make a variable that has a value of 1 any time S is below S_crit, and is zero otherwise
down1BL <- rep(0,length(s2BL)) 
down2BL <- rep(0,length(s2BL)) 
down1BL[IIBL] <- 1
down2BL[HHBL] <- 1


#make a variable that has a value of 1 any time S crosses below S_crit, and is zero otherwise
crossing1BL <- rep(0,length(s2BL)-1)
crossing2BL <- rep(0,length(s2BL)-1)
for (i in 2:length(s2BL)) {
  #if S was above S_crit in the previous timestep, but is below S_crit in the current
  #timestep, a downcrossing occurred - set "crossing" = 1 for this timestep
  if (down1BL[i-1] ==0 & down1BL[i] ==1) {  
    crossing1BL[i] <- 1    #downcrossing
  }
  if (down2BL[i-1] ==0 & down2BL[i] ==1) {  
    crossing2BL[i] <- 1    #downcrossing
  }
}

num_down1BL <- sum(crossing1BL)    #total number of downcrossings
num_down2BL <- sum(crossing2BL) 

#mean length of downcrossings (units: days)
T_bar1BL <- sum(down1BL) / num_down1BL
T_bar2BL <- sum(down2BL) / num_down2BL
#number of days that S was below S_crit divided by number of downcrossing



#number of downcrossings
IITS <- which(s2TS < swS)    #find which timesteps have an S value below S_crit
HHTS <- which(s2TS < saS)

#make a variable that has a value of 1 any time S is below S_crit, and is zero otherwise
down1TS <- rep(0,length(s2TS)) 
down2TS <- rep(0,length(s2TS)) 
down1TS[IITS] <- 1
down2TS[HHTS] <- 1


#make a variable that has a value of 1 any time S crosses below S_crit, and is zero otherwise
crossing1TS <- rep(0,length(s2TS)-1)
crossing2TS <- rep(0,length(s2TS)-1)
for (i in 2:length(s2TS)) {
  #if S was above S_crit in the previous timestep, but is below S_crit in the current
  #timestep, a downcrossing occurred - set "crossing" = 1 for this timestep
  if (down1TS[i-1] ==0 & down1TS[i] ==1) {  
    crossing1TS[i] <- 1    #downcrossing
  }
  if (down2TS[i-1] ==0 & down2TS[i] ==1) {  
    crossing2TS[i] <- 1    #downcrossing
  }
}

num_down1TS <- sum(crossing1TS)    #total number of downcrossings
num_down2TS <- sum(crossing2TS) 

#mean length of downcrossings (units: days)
T_bar1TS <- sum(down1TS) / num_down1TS
T_bar2TS <- sum(down2TS) / num_down2TS
#number of days that S was below S_crit divided by number of downcrossing



#number of downcrossings
IITL <- which(s2TL < swL)    #find which timesteps have an S value below S_crit
HHTL <- which(s2TL < saL)

#make a variable that has a value of 1 any time S is below S_crit, and is zero otherwise
down1TL <- rep(0,length(s2TL)) 
down2TL <- rep(0,length(s2TL)) 
down1TL[IITL] <- 1
down2TL[HHTL] <- 1


#make a variable that has a value of 1 any time S crosses below S_crit, and is zero otherwise
crossing1TL <- rep(0,length(s2TL)-1)
crossing2TL <- rep(0,length(s2TL)-1)
for (i in 2:length(s2TL)) {
  #if S was above S_crit in the previous timestep, but is below S_crit in the current
  #timestep, a downcrossing occurred - set "crossing" = 1 for this timestep
  if (down1TL[i-1] ==0 & down1TL[i] ==1) {  
    crossing1TL[i] <- 1    #downcrossing
  }
  if (down2TL[i-1] ==0 & down2TL[i] ==1) {  
    crossing2TL[i] <- 1    #downcrossing
  }
}

num_down1TL <- sum(crossing1TL)    #total number of downcrossings
num_down2TL <- sum(crossing2TL) 

#mean length of downcrossings (units: days)
T_bar1TL <- sum(down1TL) / num_down1TL
T_bar2TL <- sum(down2TL) / num_down2TL
#number of days that S was below S_crit divided by number of downcrossing

num_down1BS
num_down2BS
T_bar1BS
T_bar2BS

num_down1BL
num_down2BL
T_bar1BL
T_bar2BL

num_down1TS
num_down2TS
T_bar1TS
T_bar2TS

num_down1TL
num_down2TL
T_bar1TL
T_bar2TL

T_bar2BL=0

bar1 <- matrix(c(num_down1BS, num_down1BL, num_down1TS, num_down1TL, num_down2BS, num_down2BL, num_down2TS, num_down2TL), nrow=2, byrow=TRUE)      
colnames(bar1) <- c("BS","BL","TS", "TL")
rownames(bar1) <- c("num_down1", "num_down2")
bar <- as.table(bar1)
barplot(bar1,legend=T,beside=T,main='Number of Downcrossings', col=c("red", "coral"))


#bar <- matrix(c(num_down1BS, num_down1BL, num_down1TS, num_down1TL, num_down2BS, num_down2BL, num_down2TS, num_down2TL, T_bar1BS, T_bar1BL, T_bar1TS, T_bar1TL, T_bar2BS, T_bar2BL, T_bar2TS, T_bar2TL), nrow=4, byrow=TRUE)
#rownames(bar) <- c("num_down1", "num_down2")

bar2 <- matrix(c( T_bar1BS, T_bar1BL, T_bar1TS, T_bar1TL, T_bar2BS, T_bar2BL, T_bar2TS, T_bar2TL), nrow=2, byrow=TRUE)      
colnames(bar2) <- c("BS","BL","TS", "TL")
rownames(bar2) <- c("T_bar1", "T_bar2")
bar <- as.table(bar2)
barplot(bar2,legend=T,beside=T,main='Duration of Excursion', col=c("blue", "cyan"))








