#To perform ABC for MSMC:
setwd("C:/Users/Parul Johri/Work/MSMC/ABC/stats")
library("abc")
library("spatstat", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")


#t_data <- read.table("demo_disc_5_SingExon_msmc_testset_stats_sumstats_50_200.txt", h=T)
t_data <- read.table("demo_disc_5_SingExon_msmc_testset_stats_sumstats_50.txt", h=T)
t_rec <- read.table("demo_disc_5_SingExon_msmc_recombination_stats_sumstats_50.txt", h=T)
t_50 <- read.table("../../../BgsDfeDemo/abc/demo_disc_5_SingExon_sumstats_50.txt", h=T)




#get only functional stats from main file:
t_par <- t_50[,c(6:7)]
col_func <- seq(11, 76, by=3) #only linked and only means
t_stats_func <- t_50[,col_func]
t_obs_func <- t_data[,col_func]
t_rec_func <- t_rec[,col_func]

#Inference
#For functional stats only:
abc_nnet1 <- abc(target=t_obs_func[5,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet2 <- abc(target=t_obs_func[10,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet3 <- abc(target=t_obs_func[15,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet4 <- abc(target=t_obs_func[5,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet5 <- abc(target=t_obs_func[10,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet6 <- abc(target=t_obs_func[15,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")

abc_nnet_rec1 <- abc(target=t_rec_func[1,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_rec2 <- abc(target=t_rec_func[2,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_rec3 <- abc(target=t_rec_func[3,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_rec4 <- abc(target=t_rec_func[4,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_rec5 <- abc(target=t_rec_func[5,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_rec6 <- abc(target=t_rec_func[6,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
summary(abc_nnet)
summary(abc_nnet_neu)



#Plotting
par(mfcol=c(2,3))
par(mar=c(2,3,2,1))#par(mar=c(4,4,2,1))
#Equilibrium:
Nanc<-10000
Ncur<-10000
xlab <- c(expression(italic("N")["anc"]), expression(italic("N")["cur"]), expression(italic("N")["anc"]), expression(italic("N")["cur"]))
boxplot(abc_nnet_rec1$adj.values[,1], abc_nnet_rec1$adj.values[,2], abc_nnet1$adj.values[,1], abc_nnet1$adj.values[,2], ylim=c(-50,50000), cex.axis=1.8, names=xlab, border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray", lty=2)
boxplot(abc_nnet_rec4$adj.values[,1], abc_nnet_rec4$adj.values[,2], abc_nnet4$adj.values[,1], abc_nnet4$adj.values[,2], ylim=c(-50,50000), cex.axis=1.8, names=xlab, border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray", lty=2)

#Growth:
Nanc<-10000
Ncur<-20000
boxplot(abc_nnet_rec2$adj.values[,1], abc_nnet_rec2$adj.values[,2], abc_nnet2$adj.values[,1], abc_nnet2$adj.values[,2], cex.axis=1.8, names=xlab, ylim=c(-1000,50000), border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray")
boxplot(abc_nnet_rec5$adj.values[,1], abc_nnet_rec5$adj.values[,2], abc_nnet5$adj.values[,1], abc_nnet5$adj.values[,2], cex.axis=1.8, names=xlab, ylim=c(-1000,50000), border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray")

#Decline:
Nanc<-10000
Ncur<-5000
boxplot(abc_nnet_rec3$adj.values[,1], abc_nnet_rec3$adj.values[,2], abc_nnet3$adj.values[,1], abc_nnet3$adj.values[,2], cex.axis=1.8, names=xlab, ylim=c(0,50000), border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray")
boxplot(abc_nnet_rec6$adj.values[,1], abc_nnet_rec6$adj.values[,2], abc_nnet6$adj.values[,1], abc_nnet6$adj.values[,2], cex.axis=1.8, names=xlab, ylim=c(0,50000), border=c("green", "light green", "red", "red"))
abline(h=Nanc, col="black")
abline(h=Ncur, col="gray")

