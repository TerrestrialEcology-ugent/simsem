##################

# R script to reproduce the results in the submitted manuscript:

# How robust are Structural Equation Models to model mis-specification? A simulation study

# author: Lionel Hertzog

# date: 17.10.2018

########

# load the libraries
library(piecewiseSEM) # should be a version below 2.0 !!!!!!
library(lavaan)
library(plyr)
library(dplyr)
library(ggplot2)
library(igraph)
library(dagitty)
library(AICcmodavg)

# set the wd and load the helper functions
path_to_function <- "~/Desktop/Blog/GitFolder/SimSEM/"
setwd(path_to_function)

source("simSEM_function.R")

# create the parameter set to run the simulations on
# test over different data generation type (Type), sample size (N), number of variables (X),
# varying signal strength (sd_eff) and noise importance (sd_red).

sims <- expand.grid(Type=c("random","exact","shuffled","overspecified","underspecified"),
                    N=c(seq(20,100,20),200,500,1000,5000,10000),
                    X=c(5,7,10),
                    C=0.3,
                    lv=c(FALSE,TRUE),
                    sd_eff = c(1, 2.5, 5),
                    sd_res = c(0.5, 1, 2))


# for the random data type, no need to vary sd_eff idependenlty
sims <- sims[!(sims$Type == "random" & sims$sd_eff %in% c(2.5,5)),]

# run the simulations, see the code in the simSEM_function 
res_pre <- sim_for(sims,n_rep=100)

# write the output
write.table(res_pre,"simSem_res_pre.csv",row.names=FALSE)


# aggregate over the replications per parameter set
res_pre %>%
  group_by(pkg,type,Exp,N,X,Sd_eff,Sd_res)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=median(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE),
            AIC=median(AIC,na.rm=TRUE),BIC=median(BIC,na.rm=TRUE),
            hbic = median(hbic,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) -> res_all

# log-modulus transformation for thr IC metrics
res_all$AIC <- sign(res_all$AIC) * log(abs(res_all$AIC) + 1)
res_all$BIC <- sign(res_all$BIC) * log(abs(res_all$BIC) + 1)
res_all$hbic <- sign(res_all$hbic) * log(abs(res_all$hbic) + 1)

# the main figures:

## fig.1 variation in proportion of accepted models
gg1 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropSigM,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="Proportion of accepted model for different data generation")

ggsave("prop_acc.png",gg1,width=30,height=15,unit="cm",dpi=100)

## fig.2 variation in proportion of accepted paths
gg2 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropSigPaths,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  geom_hline(yintercept = 0.8,linetype="dashed") +
  labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
       title="Proportion of significant paths for different data generation")

ggsave("prop_path.png",gg2,width=30,height=15,unit="cm",dpi=100)

## fig.3 variation in average R-square
gg3 <- ggplot(subset(res_all,Sd_res == 1 & Sd_eff == 2.5),aes(x=N,y=AvgRSq,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Average R-square of the constituing models",
       title="Average R-square for different data generation")
#
ggsave("avg_rsq.png",gg3,width=30,height=15,unit="cm",dpi=100)

## fig.4 variation in proportion of failed conditional independence tests
gg4 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropCond,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Average proportion of breached conditional independence",
       title="Average failed conditional independence for different data generation")

ggsave("prop_cond.png",gg4,width=30,height=15,unit="cm",dpi=100)

## AIC figures, not shown in the main text
gg5 <- ggplot(res_all,aes(x=N,y=AIC,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Median AIC values",
       title="AIC values across simulation types")

ggsave("AIC_A.eps",gg5,width=30,height=15,unit="cm",dpi=100)

## BIC figure not shown in the main text
gg6 <- ggplot(res_all,aes(x=N,y=BIC,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Median BIC values",
       title="BIC values across simulation types")

ggsave("BIC_A.eps",gg6,width=30,height=15,unit="cm",dpi=100)

## fig.5 variation in hbic
gg7 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=hbic,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Median hbic values",
       title="HBIC values across simulation types")

ggsave("HBIC.png",gg7,width=30,height=15,unit="cm",dpi=100)

# appendix figures
## this time looking across signal / noise ratios

gg1 <- ggplot(res_all,aes(x=N,y=PropSigM,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="Proportion of accepted model for different data generation")

ggsave("prop_acc_A.eps",gg1,width=30,height=15,unit="cm",dpi=100)

gg2 <- ggplot(res_all,aes(x=N,y=PropSigPaths,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
       title="Proportion of significant paths for different data generation")

ggsave("prop_path_A.eps",gg2,width=30,height=15,unit="cm",dpi=100)
#
gg3 <- ggplot(subset(res_all,type!="random"),aes(x=N,y=AvgRSq,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Average R-square of the constituing models",
       title="Average R-square for different data generation")
#
ggsave("avg_rsq_A.eps",gg3,width=30,height=15,unit="cm",dpi=100)
#
gg4 <- ggplot(res_all,aes(x=N,y=PropCond,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Average proportion of breached conditional independence",
       title="Average failed conditional independence for different data generation")

ggsave("prop_cond_A.eps",gg4,width=30,height=15,unit="cm",dpi=100)

gg5 <- ggplot(res_all,aes(x=N,y=AIC,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Median AIC values",
       title="AIC values across simulation types")

ggsave("AIC_A.eps",gg5,width=30,height=15,unit="cm",dpi=100)

gg6 <- ggplot(res_all,aes(x=N,y=BIC,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Median BIC values",
       title="BIC values across simulation types")

ggsave("BIC_A.eps",gg6,width=30,height=15,unit="cm",dpi=100)

gg7 <- ggplot(res_all,aes(x=N,y=hbic,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Median hbic values",
       title="HBIC values across simulation types")

ggsave("HBIC_A.eps",gg6,width=30,height=15,unit="cm",dpi=100)



# compute the number of time that each scenario leads to the lowest AICc score with a difference of at least 2 to the next smallest value
lims <- c(seq(1,300,5),seq(301,780,4),seq(781,1080,5),seq(1081,1560,4),seq(1561,1860,5),seq(1861,2340,4),2341)
out <- NULL
for(i in 1:540){
  tmp <- subset(res_pre,Exp %in% (lims[i]:(lims[i+1] - 1)))
  tmp$id <- rep(1:100,nrow(tmp)/100)
  tmp_d <- dcast(tmp,id~type,value.var="AIC")
  tmp_d <- tmp_d[rowSums(tmp_d[,-1],na.rm=TRUE) != 0,] #to remove rows with all NAs
  best <- apply(tmp_d[,-1],1,function(x){
    diff <- sort(x)[2] - sort(x)[1]
    val <- ifelse(diff >= 2,
                  names(which.min(x)),
                  "unassigned")
    return(val)
  })
  preo <- as.data.frame(table(best))
  preo <- cbind(preo,tmp[1,c("N","X","Sd_eff","Sd_res","pkg")])
  names(preo)[1] <- "Type"
  preo$Exp <- i
  out <- rbind(out,preo)
  print(i)
}

#now per N-X-Sd_eff-Sd_res combination get the proportion of times each scenario got the lowest AIC score
out %>%
  group_by(Exp) %>%
  mutate(Prop = Freq / sum(Freq),Signal = paste("Sd_eff:",Sd_eff,"\nSd_res:",Sd_res)) -> out_dd

ggplot(out_dd,aes(x=N,y=Prop,color=Type,linetype=pkg)) +
  geom_path() +
  facet_grid(X~Signal) +
  scale_x_log10()

gg8 <- ggplot(subset(out_dd,Sd_eff==2.5 & Sd_res == 1),aes(x=N,y=Prop,color=Type)) +
  geom_path() +
  facet_grid(X~pkg) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of simulation that chosed a specific scenario based on AIC",
       title="AIC-based model selection")

#now the same for BIC
out <- NULL
for(i in 1:540){
  tmp <- subset(res_pre,Exp %in% (lims[i]:(lims[i+1] - 1)))
  tmp$id <- rep(1:100,nrow(tmp)/100)
  tmp_d <- dcast(tmp,id~type,value.var="BIC")
  tmp_d <- tmp_d[rowSums(tmp_d[,-1],na.rm=TRUE) != 0,] #to remove rows with all NAs
  best <- apply(tmp_d[,-1],1,function(x){
    diff <- sort(x)[2] - sort(x)[1]
    val <- ifelse(diff >= 2,
                  names(which.min(x)),
                  "unassigned")
    return(val)
  })
  preo <- as.data.frame(table(best))
  preo <- cbind(preo,tmp[1,c("N","X","Sd_eff","Sd_res","pkg")])
  names(preo)[1] <- "Type"
  preo$Exp <- i
  out <- rbind(out,preo)
  print(i)
}

#now per N-X-Sd_eff-Sd_res combination get the proportion of times each scenario got the lowest AIC score
out %>%
  group_by(Exp) %>%
  mutate(Prop = Freq / sum(Freq),Signal = paste("Sd_eff:",Sd_eff,"\nSd_res:",Sd_res)) -> out_dd

gg9 <- ggplot(subset(out_dd,Sd_eff==2.5 & Sd_res == 1),aes(x=N,y=Prop,color=Type)) +
  geom_path() +
  facet_grid(X~pkg) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of simulation that chosed a specific scenario based on BIC",
       title="BIC-based model selection")

#now for HBIC
out <- NULL
for(i in 1:540){
  tmp <- subset(res_pre,Exp %in% (lims[i]:(lims[i+1] - 1)))
  tmp$id <- rep(1:100,nrow(tmp)/100)
  tmp_d <- dcast(tmp,id~type,value.var="hbic")
  tmp_d <- tmp_d[rowSums(tmp_d[,-1],na.rm=TRUE) != 0,] #to remove rows with all NAs
  best <- apply(tmp_d[,-1],1,function(x){
    diff <- sort(abs(x))[2] - sort(abs(x))[1]
    val <- ifelse(diff >= 2,
                  names(which.min(abs(x))),
                  "unassigned")
    return(val)
  })
  preo <- as.data.frame(table(best))
  preo <- cbind(preo,tmp[1,c("N","X","Sd_eff","Sd_res","pkg")])
  names(preo)[1] <- "Type"
  preo$Exp <- i
  out <- rbind(out,preo)
  print(i)
}


#now per N-X-Sd_eff-Sd_res combination get the proportion of times each scenario got the lowest AIC score
out %>%
  group_by(Exp) %>%
  mutate(Prop = Freq / sum(Freq),Signal = paste("Sd_eff:",Sd_eff,"\nSd_res:",Sd_res)) -> out_dd

gg10 <- ggplot(out_dd,aes(x=N,y=Prop,color=Type,linetype=pkg)) +
  geom_path() +
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of simulations that chose a specific scenario based on HBIC",
       title="HBIC-based model selection")

ggsave("hbic_selection_A.eps",gg10,width=30,height=20,units = "cm")

gg10 <- ggplot(subset(out_dd,Sd_eff==2.5 & Sd_res ==1),aes(x=N,y=Prop,color=Type)) +
  geom_path() +
  facet_grid(X~pkg) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of simulations that chose a specific scenario based on HBIC",
       title="HBIC-based model selection")

ggsave("hbic_selection.png",gg10,width=15,height=20,units = "cm")