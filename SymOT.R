# The code below generates Figure 5 in the supplementary file, 
#    based on outputs from OTsymm_simulations.ipynb.

Epan = read.csv("/Users/huangzhen/desktop/HPC/Python/Epanechnikov.csv", header = F)
head(Epan)
n = 25 * (0:20) + 100
library(ggplot2)
# A colorblind friendly palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2", "#F0E442", "#D55E00", "#CC79A7")
df_Epan = data.frame(n = rep(n,5),
                     Power = c(Epan[,1],Epan[,2],Epan[,3],Epan[,4],Epan[,5]),
                     Methods = c(rep("CU(100%)",length(n)),rep("SU(100%)",length(n)),rep("t2 (86.4%)",length(n)),rep("CG(86.4%)",length(n)),rep("SG(86.4%)",length(n))))
f_Epan = ggplot(data=df_Epan[df_Epan$n>=300,], aes(x=n, y=Power, colour=Methods)) +
  geom_line() +
  scale_colour_manual(values=cbp1[c(1:4,6)]) +
  theme(legend.position=c(0.87,0.25),legend.key = element_blank(),legend.title=element_blank(),legend.background=element_blank()) 

f_Epan

Gaussian = read.csv("/Users/huangzhen/desktop/HPC/Python/Gaussian.csv",header = F)
df_Gaussian = data.frame(n = rep(n,6),
                     Power = c(Gaussian[,1],Gaussian[,2],Gaussian[,3],Gaussian[,4],Gaussian[,5],Gaussian[,6]),
                     Methods = c(rep("CentralUnif",length(n)),rep("SignUnif",length(n)),rep("t2",length(n)),rep("CentralGaussian",length(n)),rep("SignGaussian",length(n)),rep("SphericalGaussian",length(n))))
f_Gaussian = ggplot(data=df_Gaussian, aes(x=n, y=Power, colour=Methods)) +
  geom_line() +
  scale_colour_manual(values=cbp1) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 1))
f_Gaussian
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(f_Gaussian)

df_Gaussian = data.frame(n = rep(n,6),
                         Power = c(Gaussian[,1],Gaussian[,2],Gaussian[,3],Gaussian[,4],Gaussian[,5],Gaussian[,6]),
                         Methods = c(rep("CU(100%)",length(n)),rep("SU(100%)",length(n)),rep("t2 (95%)",length(n)),rep("CG(95%)",length(n)),rep("SG(95%)",length(n)),rep("SpG(95%)",length(n))))
f_Gaussian = ggplot(data=df_Gaussian, aes(x=n, y=Power, colour=Methods)) +
  geom_line() +
  scale_colour_manual(values=cbp1) +
  theme(legend.position=c(0.87,0.25),legend.key = element_blank(),legend.title=element_blank(),legend.background=element_blank()) 



gridExtra::grid.arrange(gridExtra::arrangeGrob(f_Epan + theme(legend.position=c(0.81,0.28)) + ggtitle("Epanechnikov distribution"),
                                               f_Gaussian + theme(legend.position=c(0.81,0.3)) + ggtitle("Gaussian distribution"),
                                               nrow=1),
                        mylegend, nrow=2,heights=c(10, 1))