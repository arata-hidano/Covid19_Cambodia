# Making preliminary figures as of 3 April 2020
library(ggplot2)
library(patchwork)
# Number of Clinical
#rm(gPlot1,gPlot2)
gc()
agegp = 12
X=epi_doNothingDurInf3[[1]][[1]]$time
gPlot1 <- ggplot()+ggtitle(paste0("(1) Number of clinical cases for age [",(agegp-1)*5,',',agegp*5,')',"\n",
                                  "(theoretical no intervention)"))+
  labs(x="Time(days)", y="Number of Ic")+theme_bw()+
  ylim(0,30000)

gPlot2 <- ggplot()+ggtitle(paste0("(2) Number of clinical cases for age [",(agegp-1)*5,',',agegp*5,')',"\n",
                                  "(baseline)"))+
  labs(x="Time(days)", y="Number of Ic")+theme_bw()+  ylim(0,30000)

color_v <- c("black","red","steelblue")
for(i in 1:nsim)
{
  for(j in 1:3){
    Y <- epi_doNothingDurInf3[[i]][[j]]$Ic[,agegp] 
    Y_base <-  epi_baseDurInf3[[i]][[j]]$Ic[,agegp]
    # Y_march <- epi_marchDurInf3[[i]]$ICU[,agegp]
    # Y_april <- epi_aprilDurInf3[[i]]$ICU[,agegp]
    gPlot1<-gPlot1+geom_line(aes_string(x=X, y=Y),color=color_v[j]) 
    gPlot2<-gPlot2+geom_line(aes_string(x=X, y=Y_base),color=color_v[j]) 
  }

  
}
G1 <- gPlot1 + gPlot2
G1
ggsave("NUm_IC_sc1_2.png",G1)
rm(gPlot1,gPlot2,G1)
rm(X)
gc()
# Number of ICU
rm(ICU_no,ICU_base)
X1=epi_doNothingDurInf3[[1]][[1]]$time
#rm(gPlot3,gPlot4)
gPlot3 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+ggtitle("(3) Number of ICU required for all age \n (theoretical no intervention)")+theme_bw()+
  ylim(0,80000)

gPlot4 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+
  ggtitle("(4) Number of  ICU required for all age \n (baseline)")+theme_bw()+
  ylim(0,80000)

gPlot5 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+
  ggtitle("(3) Number of  ICU required for all age \n PP vs Rural (no intervention)")+theme_bw()+
  ylim(0,80000)

gPlot6 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+
  ggtitle("(4) Number of  ICU required for all age \n PP vs Urban (no intervention)")+theme_bw()+
  ylim(0,80000)

gPlot7 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+
  ggtitle("(5) Number of  ICU required for all age \n Rural vs Urban (no intervention)")+theme_bw()+
  ylim(0,80000)

gPlot8 <- ggplot()+
  labs(x="Time(days)", y="Number of ICU")+
  ggtitle("(6) Number of  ICU required for all age \n Rural vs Urban (baseline)")+theme_bw()+
  ylim(0,80000)
color_v <- c("black","red","steelblue")

ICU_no <- ICU_base <- list()
k <- 0
for(i in 1:nsim)
{
  k <- k + 1
  for(j in 1:3){
   
    ICU_no[[k]] <- rowSums(epi_doNothingDurInf3[[i]][[j]]$ICU) 
    ICU_base[[k]] <- rowSums(epi_baseDurInf3[[i]][[j]]$ICU) 
    # Y_march <- epi_marchDurInf3[[i]]$ICU[,agegp]
    # Y_april <- epi_aprilDurInf3[[i]]$ICU[,agegp]
     
  }
  
}
for(i in 1:100)
{
  # gPlot3 <- gPlot3+geom_line(aes_string(x=X1, y=ICU_no[[1+(i-1)*3]]),color=color_v[1]) +
  #   geom_line(aes_string(x=X1, y=ICU_no[[2+(i-1)*3]]),color=color_v[2]) +
  #   geom_line(aes_string(x=X1, y=ICU_no[[3+(i-1)*3]]),color=color_v[3]) 
  # 
  # gPlot4 <- gPlot4+geom_line(aes_string(x=X1, y=ICU_base[[1+(i-1)*3]]),color=color_v[1]) +
  #   geom_line(aes_string(x=X1, y=ICU_base[[2+(i-1)*3]]),color=color_v[2]) +
  #   geom_line(aes_string(x=X1, y=ICU_base[[3+(i-1)*3]]),color=color_v[3]) 
  
  gPlot5 <- gPlot5+geom_line(aes_string(x=X1, y=ICU_no[[1+(i-1)*3]]),color=color_v[1],alpha=0.5) +
    geom_line(aes_string(x=X1, y=ICU_no[[2+(i-1)*3]]),color=color_v[2],alpha=0.5) 
  
  gPlot6 <- gPlot6+geom_line(aes_string(x=X1, y=ICU_no[[1+(i-1)*3]]),color=color_v[1],alpha=0.5) +
    geom_line(aes_string(x=X1, y=ICU_no[[3+(i-1)*3]]),color=color_v[3],alpha=0.5) 
  
  gPlot7 <- gPlot7+geom_line(aes_string(x=X1, y=ICU_no[[2+(i-1)*3]]),color=color_v[2],alpha=0.5) +
    geom_line(aes_string(x=X1, y=ICU_no[[3+(i-1)*3]]),color=color_v[3],alpha=0.5) 
  
  gPlot8 <- gPlot8+geom_line(aes_string(x=X1, y=ICU_base[[2+(i-1)*3]]),color=color_v[2],alpha=0.5) +
    geom_line(aes_string(x=X1, y=ICU_base[[3+(i-1)*3]]),color=color_v[3],alpha=0.5)
}
# G2 <- gPlot3 + gPlot4 
# ggsave("NUm_ICU_sc1_2.png",G2)
# rm(G2)
# gc()

G3 <- (gPlot5 + gPlot6)/(gPlot7 + gPlot8)
ggsave("NUm_ICU_sc1_2_separate.png",G3)
rm(G3)
gc()

