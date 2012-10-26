library(pwr)

tat_data <- read.csv("C:/Users/bakerm.MIW-VIKING/Dropbox/Tat statistics/Spots_20120920.csv", header=F)
tat_data[which(tat_data=="x")] <- NA

ind_BME <- which(tat_data[,1]=="BME")
ind_CueO <- which(tat_data[,1]=="CueO")
ind_CCCP <- which(tat_data[,1]=="CCCP")
ind_ELV <- which(tat_data[,1]=="-")

tat_data_BME <- tat_data[ind_BME,]
tat_data_CueO <- tat_data[ind_CueO,]
tat_data_CCCP <- tat_data[ind_CCCP,]
tat_data_ELV <- tat_data[ind_ELV,]

spots_CCCP <- as.numeric(as.character(tat_data_CCCP[,3]))
spots_BME <- as.numeric(as.character(tat_data_BME[,3]))
spots_CueO <- as.numeric(as.character(tat_data_CueO[,3]))
spots_ELV <- as.numeric(as.character(tat_data_ELV[,3]))


pdf("NoOfSpots_Boxplots.pdf")
boxplot(spots_CueO,spots_CCCP,spots_BME,spots_ELV, main="Distributions of Number of Spots", names=c("CueO","CCCP","BME","-"))
dev.off()

ttestRes <- t.test(spots_CueO, spots_BME)

pdf("NoOfSpots_Boxplots_CueO_BME.pdf")
boxplot(spots_CueO,spots_BME, main=paste("T-test P-value =",signif(ttestRes$p.value,3)), names=c("CueO","BME"))
dev.off()

# pdf("qqplot_CCCP.pdf")
# qqnorm(spots_CCCP)
# dev.off()
# 
# pdf("qqplot_ELV.pdf")
# qqnorm(spots_ELV)
# dev.off()

## begin reading wild-type J16 data +/- CCCP

J16_data <- read.csv("C:/Users/bakerm.MIW-VIKING//Dropbox/Tat statistics/J16_spotcounts.csv", header=F)

ind_J16 <- which(J16_data[,1]=="J16-lAry")
ind_J16_CCCP <- which(J16_data[,1]=="J16-lAry + CCCP")

J16_data_J16 <- J16_data[ind_J16,]
J16_data_CCCP <- J16_data[ind_J16_CCCP,]

spots_J16 <- as.numeric(as.character(J16_data_J16[,3]))
spots_J16_CCCP <- as.numeric(as.character(J16_data_CCCP[,3]))

t.test(spots_J16,spots_J16_CCCP)

## load in mutant data for power testing

mutant_data <- read.csv("C:/Users/bakerm.MIW-VIKING//Dropbox/Tat statistics/120924_mutant_spotcount.csv", header=F)

E16 = vector('list')
E16_F39A = vector('list')
E16_Q8A = vector('list')
E16_F39A = vector('list')
J16_F39A = vector('list')
J16_Q8A = vector('list')
J1M1 = vector('list')
J16 = vector('list')

# Put all the data into lists where list[[1]] is without CueO on plasmid, 
# list[[2]] is with CueO overexpressed and list[[3]] is with CCCP (where applicable)

E16_F39A[['-']] <- mutant_data[which(mutant_data[,1]=="E16-lAry F39A"),3]
E16_F39A[['CueO']] <- mutant_data[which(mutant_data[,1]=="E16-lAry F39A pCueO"),3]
E16_Q8A[['-']] <- mutant_data[which(mutant_data[,1]=="E16-lAry Q8A"),3]
E16_Q8A[['CueO']] <- mutant_data[which(mutant_data[,1]=="E16-lAry Q8A pCueO"),3]
J16_F39A[[1]] <- mutant_data[which(mutant_data[,1]=="J16-lAry F39A"),3]
J16_Q8A[[1]] <- mutant_data[which(mutant_data[,1]=="J16-lAry Q8A"),3]
J1M1[['-']] <- mutant_data[which(mutant_data[,1]=="J1M1-lAry"),3]
J1M1[['CueO']] <- mutant_data[which(mutant_data[,1]=="J1M1-lAry pCueO"),3]
MDABC_spots <- mutant_data[which(mutant_data[,1]=="MDABC-lAry pCueO"),3]

J16[['-']] <- spots_J16
J16[['CCCP']] <- spots_CCCP

E16[['ELV']] <- spots_ELV
E16[['CueO']] <- spots_CueO
E16[['CCCP']] <- spots_CCCP

# power testing for Q8A (eg)

power_test_spots <- function(x,y){
  
  n1 <- length(x)
  n2 <- length(y)
  mean1 <- mean(x)
  mean2 <- mean(y)
  
  pooled_sd <- sqrt( ((n1-1)*sd(x) + (n2-1)*sd(y))/(n1+n2) )
  
  d <- (mean1 - mean2) / pooled_sd
  
  pwrtest <- pwr.t.test(d = d, sig.level = 0.05, power = 0.8)
  n_needed <- pwrtest$n

}

# run power tests in various combinations. This calculates different from current N to power test
# idealised N and saves to a vector nXXX, where nXX[1] is 1 vs 2 (no CueO vs CueO)
# nXX[2] is 1 vs 3 (no CueO vs CCCP), and nXX[3] is 2 vs 3 (CueO vs CCCP)

# E16
nE16 = rep(NA,3)
pE16 = rep(NA,3)
nE16['ELV vs CueO'] <- min(length(E16[['ELV']]),length(E16[['CueO']])) - power_test_spots(E16[['ELV']],E16[['CueO']])
nE16['ELV vs CCCP'] <- min(length(E16[['ELV']]),length(E16[['CCCP']])) - power_test_spots(E16[['ELV']],E16[['CCCP']])
nE16['CueO vs CCCP'] <- min(length(E16[['CueO']]),length(E16[['CCCP']])) - power_test_spots(E16[['CueO']],E16[['CCCP']])

ttestRes <- t.test(E16[['ELV']], E16[['CueO']])
pE16['ELV vs CueO'] = ttestRes$p.value
ttestRes <- t.test(E16[['ELV']], E16[['CCCP']])
pE16['ELV vs CCCP'] = ttestRes$p.value
ttestRes <- t.test(E16[['CueO']], E16[['CCCP']])
pE16['CueO vs CCCP'] = ttestRes$p.value

# E16_F39A
nE16_F39A = vector('list')
pE16_F39A = vector('list')
nE16_F39A['- vs CueO'] <- min(length(E16_F39A[['-']]),length(E16_F39A[['CueO']])) - power_test_spots(E16_F39A[['-']],E16_F39A[['CueO']])
nE16_F39A['- vs ELV'] <- min(length(E16_F39A[['-']]),length(E16[['ELV']])) - power_test_spots(E16_F39A[['-']],E16[['ELV']])
nE16_F39A['F39A CueO vs CueO'] <- min(length(E16_F39A[['CueO']]),length(E16[['CueO']])) - power_test_spots(E16_F39A[['CueO']],E16[['CueO']])

ttestRes <- t.test(E16_F39A[['-']], E16_F39A[['CueO']])
pE16_F39A['- vs CueO'] = ttestRes$p.value
ttestRes <- t.test(E16_F39A[['-']], E16[['ELV']])
pE16_F39A['- vs ELV'] = ttestRes$p.value
ttestRes <- t.test(E16_F39A[['CueO']], E16[['CueO']])
pE16_F39A['F39A CueO vs CueO'] = ttestRes$p.value

# E16_Q8A
nE16_Q8A = vector('list')
pE16_Q8A = vector('list')
nE16_Q8A['- vs CueO'] <- min(length(E16_Q8A[['-']]),length(E16_Q8A[['CueO']])) - power_test_spots(E16_Q8A[['-']],E16_Q8A[['CueO']])
nE16_Q8A['- vs ELV'] <- min(length(E16_Q8A[['-']]),length(E16[['ELV']])) - power_test_spots(E16_Q8A[['-']],E16[['ELV']])
nE16_Q8A['Q8A CueO vs CueO'] <- min(length(E16_Q8A[['CueO']]),length(E16[['CueO']])) - power_test_spots(E16_Q8A[['CueO']],E16[['CueO']])

ttestRes <- t.test(E16_Q8A[['-']], E16_Q8A[['CueO']])
pE16_Q8A['- vs CueO'] = ttestRes$p.value
ttestRes <- t.test(E16_Q8A[['-']], E16[['ELV']])
pE16_Q8A['- vs ELV'] = ttestRes$p.value
ttestRes <- t.test(E16_Q8A[['CueO']], E16[['CueO']])
pE16_Q8A['Q8A CueO vs CueO'] = ttestRes$p.value
# J16_F39A and J16_Q8A
nJ16_mut = vector('list')
pJ16_mut = vector('list')
nJ16_mut['F39A vs J16'] <- min(length(J16_F39A[[1]]),length(J16[['-']])) - power_test_spots(J16_F39A[[1]],J16[['-']])
nJ16_mut['Q8A vs J16'] <- min(length(J16_Q8A[[1]]),length(J16[['-']])) - power_test_spots(J16_Q8A[[1]],J16[['-']])


ttestRes <- t.test(J16_F39A[[1]], J16[['-']])
pJ16_mut['F39A vs J16'] = ttestRes$p.value
ttestRes <- t.test(J16_Q8A[[1]], J16[['-']])
pJ16_mut['Q8A vs J16'] = ttestRes$p.value

#J1M1
nJ1M1 = vector('list')
pJ1M1 = vector('list')
nJ1M1['- vs ELV'] <- min(length(J1M1[['-']]),length(E16[['ELV']])) - power_test_spots(J1M1[['-']],E16[['ELV']])
nJ1M1['CueO vs CueO'] <- min(length(J1M1[['CueO']]),length(E16[['CueO']])) - power_test_spots(J1M1[['CueO']],E16[['CueO']])
nJ1M1['- vs CueO'] <- min(length(J1M1[['-']]),length(J1M1[['CueO']])) - power_test_spots(J1M1[['-']],J1M1[['CueO']])

ttestRes <- t.test(J1M1[['-']], E16[['ELV']])
pJ1M1['- vs ELV'] = ttestRes$p.value
ttestRes <- t.test(J1M1[['CueO']], E16[['CueO']])
pJ1M1['CueO vs CueO'] = ttestRes$p.value
ttestRes <- t.test(J1M1[['-']], J1M1[['CueO']])
pJ1M1['- vs CueO'] = ttestRes$p.value
#MDABC


####
#mutants from 20121009

tat_data <- read.csv("C:/Users/bakerm.MIW-VIKING/Dropbox/Tat statistics/spotcount_20121009-HiPip+DADE.csv", header=F)

DADE_F39A = vector('list')
E16_HP = vector('list')

DADE_F39A[['-']] <- tat_data[which(tat_data[,1]=="DADE-lAry F39A"),3]
E16[['ELV 2']] <- tat_data[which(tat_data[,1]=="E16-lAry"),3]
E16[['CCCP 2']] <- tat_data[which(tat_data[,1]=="E16-lAry + CCCP"),3]
E16_F39A[['-2']] <- tat_data[which(tat_data[,1]=="E16-lAry F39A"),3]
E16_F39A[['CCCP2']] <- tat_data[which(tat_data[,1]=="E16-lAry F39A + CCCP"),3]
E16_HP[['HiPip']] <- tat_data[which(tat_data[,1]=="ElAry pHP"),3]
E16_HP[['CCCP']] <- tat_data[which(tat_data[,1]=="ElAry pHP CCCP"),3]
E16_HP[['bME']] <- tat_data[which(tat_data[,1]=="ElAry pHP CCCP bME"),3]
J1M1[['p101BC']] <- tat_data[which(tat_data[,1]=="J1M1 p101BC pCueO"),3]
J1M1[['p101*BC']] <- tat_data[which(tat_data[,1]=="J1M1 p101*BC pCueO"),3]
J1M1[['F39A']] <- tat_data[which(tat_data[,1]=="J1M1-lAry F39A"),3]

hist
##save boxplot pdfs

pdf("BoxPlot_J1M1.pdf")
boxplot(J1M1, main=paste("J1M1"))
dev.off()
pdf("BoxPlot_E16.pdf")
boxplot(E16, main=paste("E16"))
dev.off()
pdf("BoxPlot_E16_F39A.pdf")
boxplot(E16_F39A, main=paste("E16_F39A"))
dev.off()
pdf("BoxPlot_E16_Q8A.pdf")
boxplot(E16_Q8A, main=paste("E16_Q8A"))
dev.off()
pdf("BoxPlot_E16HP.pdf")
ttestRes <- t.test(E16_HP[['HiPip']],E16_HP[['bME']])
boxplot(E16_HP, main=paste("T-test P-value =",signif(ttestRes$p.value,3)))
dev.off()

## putting n calculations and p values inside lists

nE16_F39A['- vs CCCP 2'] <- min(length(E16_F39A[['-2']]),length(E16_F39A[['CCCP2']])) - power_test_spots(E16_F39A[['-2']],E16_F39A[['CCCP2']])
nE16['ELV2 vs CCCP2'] <- min(length(E16[['ELV 2']]),length(E16[['CCCP 2']])) - power_test_spots(E16[['ELV 2']],E16[['CCCP 2']])
nJ1M1['CueO vs p101*BC'] <- min(length(J1M1[['CueO']]),length(J1M1[['p101*BC']])) - power_test_spots(J1M1[['CueO']],J1M1[['p101*BC']])
nJ1M1['p101BC vs p101*BC'] <- min(length(J1M1[['p101BC']]),length(J1M1[['p101*BC']])) - power_test_spots(J1M1[['p101BC']],J1M1[['p101*BC']])
nE16['pCueO vs F39A'] <- min(length(E16[['CueO']]),length(E16_F39A[['-']])) - power_test_spots(E16[['CueO']],E16_F39A[['-']])
pE16['pCueO vs F39A'] <-  t.test(E16[['CueO']],E16_F39A[['-']])$p.value
nJ1M1['pCueO vs F39A']<- min(length(J1M1[['CueO']]),length(J1M1[['F39A']])) - power_test_spots(J1M1[['CueO']],J1M1[['F39A']])
pJ1M1['pCueO vs F39A'] <- t.test(J1M1[['CueO']],J1M1[['F39A']])$p.value

rough_poisson_n <- function(X1,X2) {
  
  n <- 4/((sqrt(mean(X1)) - sqrt(mean(X2)))^2)

}

# wilcox.test(J1M1[['CueO']],J1M1[['F39A']])