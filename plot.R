library("scales")
library("ggthemes")
library("ggsci")

grouptest<-survfit(Surv(time,delta)~r_group_ind, data=groupdata)

par(mar=c(5.1, 4.1, 4.1, 3.1))
summ<-summary(grouptest,time=60)
plot(grouptest,ylim=c(0, 1.0),
     
     xlim=c(0,60),
     xlab=("Survival Time in Months"),
     ylab=("Proportion Surviving"),
     cex=1.2,
     cex.lab=1,
     cex.axis=1,
     col=gdocs_pal()(10)[2:9],
     lwd=2.0,
     lty=c(2,2,2,2,1,1,1,1),
     mark.time=FALSE)
title(main = "")
legend("bottomleft",
       legend=legend,
       col=gdocs_pal()(10)[2:9],
       cex=0.7,
       lty=c(2,2,2,2,1,1,1,1),
       lwd=2.0,
       xjust=1, yjust=1)


rysurv<-paste(round(100*summ$surv, 1), "%", sep="")
text(66,cex = 0.6,summ$surv,rysurv,col=gdocs_pal()(10)[2:9], xpd=NA)

par(mar=c(5.1, 4.1, 4.1, 2.1))


##################################### Dendrogram ############################################
library(factoextra)

# Plot C-indices against number of groups and give n*
weights1<-1:numComb
for (i in 1:length(assign)){
  weights1[unlist(assign[i])]<-100^i
}

weights2<-rank(-SP)


myClust$labels <-as.character(label)
#fviz_dend(reorder(as.dendrogram(myClust),numComb+1-rank(factor(label))),
fviz_dend(reorder(as.dendrogram(myClust),weights1+weights2),          
          repel = F,
          ylab="Dissimilarity",
          main='',
          k=finalGroups,
          lwd = 0.3,
          cex = 0.3,                     # Label size
          palette = gdocs_pal()(10)[2:9],                # Color palette see ?ggpubr::ggpar
          #k_colors = c("salmon2", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#8F3931FF", "#58593FFF", "#350E20FF"),
          rect = TRUE, 
          rect_fill = TRUE, # Add rectangle around groups
          rect_border = gdocs_pal()(10)[2:9], 
          labels_track_height = 0.5      # Augment the room for labels
)


# Plot C-indices against number of groups and give n*
plot(0,0,xlim = c(1,numComb),ylim = c(0.5,min(1,max(allC)+0.05)),cex.axis=1,type = "n",main="",xlab="Number of Groups",ylab="C-index",cex.lab=1)
lines(1:numComb, allC  , col =" darkblue ", lty =1, lwd =2)
abline(v=bestn,col='gray', lwd =1)
mtext(paste("n* =",bestn), cex = 1,adj=0.5, padj=0.2,side=1,at=bestn, col = "darkblue")
points(bestn, allC[bestn],type="p", pch=19, col="darkblue",cex = 1.4)
text(bestn, y = allC[bestn], labels =paste("EACCD: C-index = ",formatC(allC[bestn], format='f',digit=4),sep=''), adj = c(-0.05,+1.2),
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1, col = "darkblue", font = NULL)




points(24, 0.6207,type="p", pch=19, col="black",cex = 1.4)
text(24, y = 0.6207, labels =paste("Ann Arbor: C-index = ",0.6207,sep=''), adj = c(-0.05,+1.2),
     pos = NULL, offset = 2, vfont = NULL,
     cex = 1, col = "black", font = NULL)


plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
