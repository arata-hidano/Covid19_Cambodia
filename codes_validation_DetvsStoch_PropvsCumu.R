# Code validation for (1) stochastic model and (2) differnece equation for resource model
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

# Compare Deterministic and stochastic model
# ICU
sum_ICU_det <- rowSums(epi_doNothingDurInf3[[1]]$ICU_det)
X=epi_doNothingDurInf3[[1]]$time

ICU_comparison <- ggplot()+geom_line(aes_string(x=X, y=sum_ICU_det),color="black",size=1.5)+
  ggtitle("(A) Number of required ICU over time")+
  labs(x="Time(days)", y="Number of ICU")+theme_bw()


for(i in 1:nsim)
{
  sum_ICU_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$ICU)
  ICU_comparison <- ICU_comparison + geom_line(aes_string(x=X, y=sum_ICU_stochastic),color="steelblue",size=0.5,alpha=0.5)
  
}
ICU_comparison
#ggsave(file="ICU_comparison.png",ICU_comparison)
#BED
sum_BED_det <- rowSums(epi_doNothingDurInf3[[1]]$BED_det)
X=epi_doNothingDurInf3[[1]]$time

BED_comparison <- ggplot()+geom_line(aes_string(x=X, y=sum_BED_det),color="black",size=1.5)+
  ggtitle("(B) Number of required BED over time")+
  labs(x="Time(days)", y="Number of BED")+theme_bw()


for(i in 1:nsim)
{
  sum_BED_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$BED)
  BED_comparison <- BED_comparison + geom_line(aes_string(x=X, y=sum_BED_stochastic),color="steelblue",size=0.5,alpha=0.5)
  
}
BED_comparison
#ggsave(file="BED_comparison.png",BED_comparison)

# Check cumulative number of FA, BED, ICU against that of Ic
# Retrieve parameters
durDelay = 7;                                                 # Median time from symptom onset to hospital admission Wang et al 2020
delta = 0.2;                                                       # proportion of clinical cases that develop severe signs
theta = 1-exp(-1/durDelay);                                       # Probability of hospital admission after symptom onset
epsilon = 0.26 
# FA
X=epi_doNothingDurInf3[[1]]$time

cum_Ic_FA <- ggplot()+geom_hline(aes(yintercept = delta),color="black",size=0.5)+
  ggtitle("(C) Propotion of severe cases \n among clinical cases over time")+
  labs(x="Time(days)", y="Proportion")+theme_bw()


for(i in 1:nsim)
{
  sum_Ic_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$cum_Ic)
  sum_FA_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$cum_FA)
  prop <- sum_FA_stochastic/sum_Ic_stochastic
  cum_Ic_FA <- cum_Ic_FA + geom_line(aes_string(x=X, y=prop),color="steelblue",size=0.5)
  
}
cum_Ic_FA

# ICU
X=epi_doNothingDurInf3[[1]]$time

cum_ICU_FA <- ggplot()+geom_hline(aes(yintercept = epsilon),color="black",size=0.5)+
  ggtitle("(D) Propotion of ICU cases \n among severe cases over time")+
  labs(x="Time(days)", y="Proportion")+theme_bw()


for(i in 1:nsim)
{
  sum_ICU_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$cum_ICU)
  sum_FA_stochastic <- rowSums(epi_doNothingDurInf3[[i]]$cum_FA)
  prop <- sum_ICU_stochastic/sum_FA_stochastic
  cum_ICU_FA <- cum_ICU_FA + geom_line(aes_string(x=X, y=prop),color="steelblue",size=0.5)
  
}
cum_ICU_FA

combine <- (ICU_comparison+BED_comparison)/(cum_Ic_FA+cum_ICU_FA)
ggsave("Validation_StochasticDet_PropCum.png",combine)