library(pwr)

tat_data <- read.csv("C:/Users/bakerm.MIW-VIKING/Dropbox/Tat statistics/grandspotcount.csv", header=F)
tat_data[which(tat_data=="x")] <- NA

E = vector('list')
F39A = vector('list')
Q8A = vector('list')
J = vector('list')
D = vector('list')
J1M1 = vector('list')
MDABC = vector('list')

E[['pCueO']] <- tat_data[which(tat_data[,1]=="ElArY pCueO"),3]
E[['pCueO CCCP']] <- tat_data[which(tat_data[,1]=="ElArY pCueO CCCP"),3]
E[['bME']] <- tat_data[which(tat_data[,1]=="ElArY pCueO CCCP BME"),3]
E[['kk']] <- tat_data[which(tat_data[,1]=="ElArY pCueOkk"),3]

E[['pHP']] <- tat_data[which(tat_data[,1]=="ElArY pHP"),3]
E[['pHP CCCP']] <- tat_data[which(tat_data[,1]=="ElArY pHP CCCP"),3]
E[['pHP bME']] <- tat_data[which(tat_data[,1]=="ElArY pHP CCCP bME"),3]
E[['pHPkk']] <- tat_data[which(tat_data[,1]=="ElArY pHPkk"),3]

E[['ArY']] <- tat_data[which(tat_data[,1]=="ElArY"),3]
E[['ArY CCCP']] <- tat_data[which(tat_data[,1]=="ElArY CCCP"),3]
E[['kk CCCP']] <- tat_data[which(tat_data[,1]=="CueOkk CCCP"),3]



F39A[['pCueO']] <- tat_data[which(tat_data[,1]=="ElArY F39A pCueO"),3]
F39A[['CCCP']] <- tat_data[which(tat_data[,1]=="ElArY F39A CCCP"),3]
F39A[['ArY']]  <- tat_data[which(tat_data[,1]=="ElArY F39A"),3]

Q8A[['pCueO']] <- tat_data[which(tat_data[,1]=="ElArY Q8A pCueO"),3]
# Q8A[['CCCP']] <- tat_data[which(tat_data[,1]=="ElArY Q8A CCCP"),3]
Q8A[['ArY']]  <- tat_data[which(tat_data[,1]=="ElArY Q8A"),3]

J1M1[['pCueO']] <- tat_data[which(tat_data[,1]=="J1M1lArY pCueO"),3]
J1M1[['ArY']] <- tat_data[which(tat_data[,1]=="J1M1lArY"),3]
J1M1[['F39A']] <- tat_data[which(tat_data[,1]=="J1M1lArY F39A"),3]
J1M1[['pBC']] <- tat_data[which(tat_data[,1]=="J1M1lArY pBC pCueO"),3]
J1M1[['pBC*']] <- tat_data[which(tat_data[,1]=="J1M1lArY pBC* pCueO"),3]

J[['ArY']] <- tat_data[which(tat_data[,1]=="JlArY"),3]
J[['CCCP']] <- tat_data[which(tat_data[,1]=="JlArY CCCP"),3]
J[['F39A']] <- tat_data[which(tat_data[,1]=="JlArY F39A"),3]
J[['Q8A']] <- tat_data[which(tat_data[,1]=="JlArY Q8A"),3]

D[['F39A']] <- tat_data[which(tat_data[,1]=="DADElArY F39A"),3]

MDABC[['CueO']] <- tat_data[which(tat_data[,1]=="MDABC-lArY pCueO"),3]
MDABC[['ArY']] <- tat_data[which(tat_data[,1]=="MDABC-lArY"),3]

## plot histograms for each

# hist(x, breaks = "Sturges",
#      freq = NULL, probability = !freq,
#      include.lowest = TRUE, right = TRUE,
#      density = NULL, angle = 45, col = NULL, border = NULL,
#      main = paste("Histogram of" , xname),
#      xlim = range(breaks), ylim = NULL,
#      xlab = xname, ylab,
#      axes = TRUE, plot = TRUE, labels = FALSE,
#      nclass = NULL, warn.unused = TRUE, ...)

## print all histograms

for (name in names(E)) {
  pdf(paste("hist_BCE", name, ".pdf", sep="_"))
  hist(E[[name]],xlab = "Number of Spots", main = paste("Histogram of AryBCE" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE, axes = F)
  axis(2)
  axis(1, at=c(0:10))
  dev.off()
}

for (name in names(J1M1)) {
  pdf(paste("hist_ABCE", name, ".pdf", sep="_"))
  hist(J1M1[[name]],xlab = "Number of Spots", main = paste("Histogram of AAryBCE" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

for (name in names(J)) {
  pdf(paste("hist_BC", name, ".pdf", sep="_"))
  hist(J[[name]],xlab = "Number of Spots", main = paste("Histogram of AryBC" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

for (name in names(F39A)) {
  pdf(paste("hist_BCE_F39A", name, ".pdf", sep="_"))
  hist(F39A[[name]],xlab = "Number of Spots", main = paste("Histogram of AryBCE F39A" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

for (name in names(Q8A)) {
  pdf(paste("hist_BCE_Q8A", name, ".pdf", sep="_"))
  hist(Q8A[[name]],xlab = "Number of Spots", main = paste("Histogram of AryBCE Q8A", name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

for (name in names(D)) {
  pdf(paste("hist_DADE_F39A", name, ".pdf", sep="_"))
  hist(D[[name]],xlab = "Number of Spots", main = paste("Histogram of Ary F39A" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

for (name in names(MDABC)) {
  pdf(paste("hist_E F39A", name, ".pdf", sep="_"))
  hist(MDABC[[name]],xlab = "Number of Spots", main = paste("Histogram of AryE" , name),
       col = 'grey',breaks = c(0:10),xlim = range(0,10),ylim = range(0,1),freq = FALSE,axes = F)
  axis(2)
  axis(1, at=c(0:10)))
  dev.off()
}

pdf("Boxplot_BCE.pdf")
boxplot(E[['pCueO']],E[['pCueO CCCP']],E[['bME']],E[['kk']],E[['pHP']],E[['pHP CCCP']], E[['pHP bME']],E[['pHPkk']],main = paste("Boxplot of AryBCE"),xlab = "Strain", ylab = "Number of Spots")
box()
axis(2)
axis(1, at=1:length(names(E)), labels=abbreviate(names(E)))
dev.off()

# Distributions to compare that are quite similar to each other (but might be different):
#   
#   ElAry vs ElAry CCCP
# ElAry vs ElAry pCueOkk vs pHPkk
# ElAry pCueO vs ElAry pHP (NB can maybe include the BME data for each strain here to increase n??)
# ElAry pCueO vs JlAry
# ElAry F39A vs ElAry F39A CCCP
# ElAry pCueO vs ElAry Q8A pCueO

t.test(E[['ArY']],E[['ArY CCCP']])$p.value
t.test(E[['pCueO']],E[['pHP']])$p.value
t.test(E[['pCueO']],J[['ArY']])$p.value
t.test(F39A[['ArY']],F39A[['CCCP']])$p.value
t.test(E[['pCueO']],Q8A[['pCueO']])$p.value

Ary_vs_CueOkk_vs_pHPkk =c( t.test(E[['ArY']],E[['kk']])$p.value,
                           t.test(E[['ArY']],E[['pHPkk']])$p.value,t.test(E[['kk']],E[['pHPkk']])$p.value)


## now let's calculate p values for every internal combination of each list
# first we have to unlist these components, so t.test will work (as only works on numerics)
# and we do a double for loop so that we do every combination, filling out an upper triangular matrix

pval = vector('list')


pval[['E']] = matrix(data = NA, nrow = length(E), ncol = length(E), byrow = FALSE,
               dimnames = list(names(E),names(E)))
for (i in 1:length(E)) for (j in i:length(E)) pval[['E']][i,j]<- t.test(unlist(E[i]),unlist(E[j]))$p.value

pval[['J']] = matrix(data = NA, nrow = length(J), ncol = length(J), byrow = FALSE,
               dimnames = list(names(J),names(J)))
for (i in 1:length(J)) for (j in i:length(J)) pval[['J']][i,j]<- t.test(unlist(J[i]),unlist(J[j]))$p.value

pval[['J1M1']] = matrix(data = NA, nrow = length(J1M1), ncol = length(J1M1), byrow = FALSE,
                  dimnames = list(names(J1M1),names(J1M1)))
for (i in 1:length(J1M1)) for (j in i:length(J1M1)) pval[['J1M1']][i,j]<- 
  t.test(unlist(J1M1[i]),unlist(J1M1[j]))$p.value

pval[['F39A']] = matrix(data = NA, nrow = length(F39A), ncol = length(F39A), byrow = FALSE,
                  dimnames = list(names(F39A),names(F39A)))
for (i in 1:length(F39A)) for (j in i:length(F39A)) pval[['F39A']][i,j]<- 
  t.test(unlist(F39A[i]),unlist(F39A[j]))$p.value

pval[['Q8A']] = matrix(data = NA, nrow = length(Q8A), ncol = length(Q8A), byrow = FALSE,
                  dimnames = list(names(Q8A),names(Q8A)))
for (i in 1:length(Q8A)) for (j in i:length(Q8A)) pval[['Q8A']][i,j]<- 
  t.test(unlist(Q8A[i]),unlist(Q8A[j]))$p.value

signif_thresh <- 0.05
pval[['E']]<signif_thresh
pval[['J']]<signif_thresh
pval[['J1M1']]<signif_thresh
pval[['F39A']]<signif_thresh
pval[['Q8A']]<signif_thresh

bonfo = vector('list')

bonfo[['E']] <- pval[['E']]<signif_thresh/(length(E)-1)
bonfo[['J']] <- pval[['J']]<signif_thresh/(length(J)-1)
bonfo[['J1M1']] <- pval[['J1M1']]<signif_thresh/(length(J1M1) - 1)
bonfo[['F39A']] <- pval[['F39A']]<signif_thresh/(length(F39A) -1)
bonfo[['Q8A']] <- pval[['Q8A']]<signif_thresh/(length(Q8A) -1)
