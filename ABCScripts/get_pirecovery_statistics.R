#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_pirecovery_statistics.R
#Rscript ./get_pirecovery_statistics.R 3 500 geneconv_humans 50000 4000
#setwd("/home/pjohri/eqm_disc_5/")

#options("winSize" = ()$win_size)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
simID <- args[1]
win_size <- as.numeric(args[2])
#step_size <- as.numeric(args[3])
folder <- args[3]
noncoding_size <- as.numeric(args[4])
coding_size <- as.numeric(args[5])
#simID
#win_size
t <- read.table(paste("/scratch/pjohri1/", folder, "/sim", simID, "_", args[2], ".stats", sep=""),h=T,fill=T)

chr_len <- noncoding_size + coding_size
num_win = chr_len/win_size

#To calculate mean and variance of all statistics per window:
for (i in 1:num_win) {
	wpii <- t$thetapi[which(t$posn==i)]/win_size
	wwi <- t$thetaw[which(t$posn==i)]/win_size
	whi <- t$thetah[which(t$posn==i)]/win_size
	whpi <- t$hprime[which(t$posn==i)]
	wtdi <- t$tajimasd[which(t$posn==i)]
	wsingi <- t$numSing[which(t$posn==i)]/win_size
	whapdivi <- t$hapdiv[which(t$posn==i)]
	wrsqi <- as.numeric(as.character(t$rsq[which(t$posn==i & t$rsq!="<NA>")]))
	wDi <- t$D[which(t$posn==i)]
	wDpri <- t$Dprime[which(t$posn==i)]
	wdivi <- t$div[which(t$posn==i)]/win_size
	v_win_m <- c(mean(wpii, na.rm=T), mean(wwi, na.rm=T), mean(whi, na.rm=T), mean(whpi, na.rm=T), mean(wtdi, na.rm=T), mean(wsingi, na.rm=T), mean(whapdivi, na.rm=T), mean(wrsqi, na.rm=T), mean(wDi, na.rm=T), mean(wDpri, na.rm=T), mean(wdivi, na.rm=T))
	v_win_sd <- c(sd(wpii, na.rm=T), sd(wwi, na.rm=T), sd(whi, na.rm=T), sd(whpi, na.rm=T), sd(wtdi, na.rm=T), sd(wsingi, na.rm=T), sd(whapdivi, na.rm=T), sd(wrsqi, na.rm=T), sd(wDi, na.rm=T), sd(wDpri, na.rm=T), sd(wdivi, na.rm=T))
	if (i==1){
		stats <- matrix(c(v_win_m, v_win_sd))
		colnames(stats) <- c(as.character(i))
		rownames(stats) <- c("thetapi_m", "thetaw_m", "thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd")
		}
	else{
		
		stats <- cbind(stats, c(v_win_m, v_win_sd))
		}
	}
colnames(stats, do.NULL = FALSE)
colnames(stats) <- c(1:num_win)

write.table(stats, file=paste("/scratch/pjohri1/", folder, "/sim", simID, "_", args[2], ".winsummary", sep=""), sep="\t")


#Look at recovery of only pi, fit it to log curve and write out those values:
#We are fitting only the left side:
num_win_nc <- noncoding_size/win_size
x <- c(1:num_win_nc)
xx <- c(1:num_win_nc)
y <- stats["thetapi_m",][num_win_nc:1]
fit_log <- lm(y~log(x))
intercept <- fit_log$coefficients[1]
slope <- fit_log$coefficients[2]
max <- slope*log(num_win_nc)+intercept
predicted_max <- slope*log(10000)+intercept
reduction <- ((max-intercept)/max)*100
predicted_reduction <- ((predicted_max-intercept)/predicted_max)*100
theory_max <- 0.02
#number of bases needed for recovery:
y50 <- 0.5*(max - intercept)+intercept
y75 <- 0.75*(max - intercept)+intercept
y90 <- 0.90*(max - intercept)+intercept
numwin50 <- exp((y50-intercept)/slope)
numwin75 <- exp((y75-intercept)/slope)
numwin90 <- exp((y90-intercept)/slope)
numbp50 <- numwin50*win_size
numbp75 <- numwin75*win_size
numbp90 <- numwin90*win_size

v_labels <- c("slope", "intercept", "max", "predicted_max", "reduction_pi", "predicted_reduction_pi", "numbp50", "numbp75", "numbp90")
v_values <- c(slope, intercept, max, predicted_max, reduction, predicted_reduction, numbp50, numbp75, numbp90)
write(v_labels, file=paste("/scratch/pjohri1/", folder, "/sim", simID, "_", args[2], ".recoverypi", sep=""), ncolumns=length(v_labels), append = F, sep = '\t')
write(v_values, file=paste("/scratch/pjohri1/", folder, "/sim", simID, "_", args[2], ".recoverypi", sep=""), ncolumns=length(v_values), append = T, sep = '\t')


print ("done")



