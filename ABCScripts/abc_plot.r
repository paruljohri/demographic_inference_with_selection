#To perform ABC for MSMC:
setwd("C:/Users/Parul Johri/Work/MSMC/ABC/stats")
library("abc")
library("spatstat", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")


#t_data <- read.table("demo_disc_5_SingExon_msmc_testset_stats_sumstats_50_200.txt", h=T)
t_data <- read.table("demo_disc_5_SingExon_msmc_testset_stats_sumstats_50.txt", h=T)
t_50 <- read.table("../../../BgsDfeDemo/abc/demo_disc_5_SingExon_sumstats_50.txt", h=T)
t_neu <- read.table("demo_neutral_SingExon_msmc_stats_sumstats_50.txt", h=T)
t_par <- t_neu[,c(6:7)]

#get linked stats (only means) from main file:
col_link <- seq(12, 44, by=3) #only linked and only means
t_stats_link <- t_neu[,col_link]
to_remove <- c("link_div_m")
t_stats <- t_stats_link[ , !(names(t_stats_link) %in% to_remove)]

#get stats from observation file:
col_link <- seq(12, 44, by=3)
t_obs_link <- t_data[,col_link]
to_remove <- c("link_div_m")
t_obs <- t_obs_link[ , !(names(t_obs_link) %in% to_remove)]


#Check with stats that only have SFS information:
col_link <- seq(12, 28, by=3) #only linked and only means
t_stats <- t_neu[,col_link]
#get stats from observation file:
col_link <- seq(12, 28, by=3)
t_obs <- t_data[,col_link]

#get only functional stats from main file:
t_par <- t_50[,c(6:7)]
t_par_neu <- t_neu[,c(6:7)]
col_func <- seq(11, 76, by=3) #only linked and only means
t_stats_func <- t_50[,col_func]
t_obs_func <- t_data[,col_func]
t_stats_neu_func <- t_neu[,col_func]

#Inference
#For functional stats only:
abc_nnet1 <- abc(target=t_obs_func[1,], param=t_par, sumstat=t_stats_func, tol=0.1, method="neuralnet")
abc_nnet_neu1 <- abc(target=t_obs_func[1,], param=t_par_neu, sumstat=t_stats_neu_func, tol=0.1, method="neuralnet")
summary(abc_nnet)
summary(abc_nnet_neu)



#Plotting
par(mfrow=c(1,5))
par(mar=c(4,4,2,1))
#Equilibrium:
xlab <- c(expression(italic("N")["anc"]), expression(italic("N")["cur"]))
boxplot(abc_nnet_neu1$adj.values, ylim=c(0,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet_neu2$adj.values, ylim=c(0,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet_neu3$adj.values, ylim=c(0,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet_neu4$adj.values, ylim=c(0,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet_neu5$adj.values, ylim=c(0,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
>>save as A4 (4 x 11.76)

boxplot(abc_nnet1$adj.values, ylim=c(5000,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet2$adj.values, ylim=c(5000,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet3$adj.values, ylim=c(5000,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet4$adj.values, ylim=c(5000,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)
boxplot(abc_nnet5$adj.values, ylim=c(5000,15000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=10000, col="purple")
abline(h=10000, col="green", lty=2)

#Growth:
Nanc<-10000
Ncur<-20000
xlab <- c(expression(italic("N")["anc"]), expression(italic("N")["cur"]))
boxplot(abc_nnet_neu6$adj.values, cex.axis=2, names=xlab, ylim=c(-1000,30000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu7$adj.values, cex.axis=2, names=xlab, ylim=c(-1000,30000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu8$adj.values, cex.axis=2, names=xlab, ylim=c(-1000,30000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu9$adj.values, cex.axis=2, names=xlab, ylim=c(-1000,30000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu10$adj.values, cex.axis=2, names=xlab, ylim=c(-1000,30000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")

boxplot(abc_nnet6$adj.values, ylim=c(5000,25000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet7$adj.values, ylim=c(5000,25000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet8$adj.values, ylim=c(5000,25000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet9$adj.values, ylim=c(5000,25000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet10$adj.values, ylim=c(5000,25000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")

#Decline:
Nanc<-10000
Ncur<-5000
xlab <- c(expression(italic("N")["anc"]), expression(italic("N")["cur"]))
boxplot(abc_nnet_neu11$adj.values, cex.axis=2, names=xlab, ylim=c(0,11000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu12$adj.values, cex.axis=2, names=xlab, ylim=c(0,11000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu13$adj.values, cex.axis=2, names=xlab, ylim=c(0,11000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu14$adj.values, cex.axis=2, names=xlab, ylim=c(0,11000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet_neu15$adj.values, cex.axis=2, names=xlab, ylim=c(0,11000), border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")

boxplot(abc_nnet11$adj.values, ylim=c(2000,13000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet12$adj.values, ylim=c(2000,13000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet13$adj.values, ylim=c(2000,13000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet14$adj.values, ylim=c(2000,13000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")
boxplot(abc_nnet15$adj.values, ylim=c(2000,13000), cex.axis=2, names=xlab, border=c("green", "purple"))
abline(h=Nanc, col="green")
abline(h=Ncur, col="purple")