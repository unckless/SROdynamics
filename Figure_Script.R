library(ggplot2); library(cowplot); library(metR); library(dplyr); library(forcats)

##Figure 2 - invasion (idealized)
df=data.frame(matrix(NA,ncol=3,nrow=0))
names(df)=c("so","u","dYinvade")

for(so in seq(0.01,1,0.01)) {
  for(u in 10^seq(-5,-1,0.05)) {
   
      
      df[length(df$dY)+1,]=c(so,u,(u-so*u)/(2*so-2*u))
    }}

df$dYinvade[which(df$dYinvade>0.5)]=NA
df$dYinvade[which(df$dYinvade<0)]=NA

dfsimple=df #Store idealized model

## Explore other parameters for invasion
df=data.frame(matrix(NA,ncol=6,nrow=0))
names(df)=c("so","u","dO","ssro","ssry","dYinvade")

dO=-1/2; ssro=0.5;
for(so in seq(0.01,1,0.01)) {
  for(u in 10^seq(-5,-1,0.05)) {
    for(ssry in c(0,0.2,0.4)) {
        if(ssro==0) ssro=so
    df[length(df$dY)+1,]=c(so,u,dO,ssro,ssry,(-so*ssry+(so-2*dO*(-1+ssro)-ssro+ssry)*u)/(2*(-1+ssry)*(so-u)))
    
  }}}

df$dYinvade[which(df$dYinvade>0.5)]=NA
df$dYinvade[which(df$dYinvade<0)]=NA

dfsimple$ssry=-1
df2=rbind(dfsimple,df[,c(1,2,6,5)])

pdf("lambda_SR_invade.pdf")
ggplot(df2,aes(x=so,y=u,z=dYinvade,fill=dYinvade))+geom_tile()+scale_y_continuous(trans='log10')+theme_cowplot()+
  geom_contour(col="darkgray",na.rm=TRUE,breaks=seq(0,0.5,0.1))+
  geom_text_contour(aes(z=dYinvade),breaks=seq(0,0.5,0.1),col="darkgray",nudge_x=0.03,nudge_y=-0.02,rotate=FALSE)+
  scale_fill_gradient(low="white",high="red",na.value="gray")+facet_wrap(.~ssry,ncol=2)
dev.off()

###Figure 3A and B
df=data.frame(matrix(NA,ncol=5,nrow=0))
names(df)=c("dY","lambda","so","Chrom","Eq")


u=0.001
for(so in c(0.0,0.2,0.4,0.8)) {
  for(dY in seq(0.05,0.5,0.01)) {
    for(lambda in 10^seq(0,3,0.1)) {
      
      df[length(df$dY)+1,]=c(dY,lambda,so,"O",1-(1/(2+8*dY^2-8*dY*(-1+so)-2*so*(2+so)))*
                               (2-2*so-2*so^2+u+so*u+lambda*u-so*lambda*u-4*dY^2*(-1+lambda*u)+2*dY*(3+u+so*(-3+lambda*u))-
                                  sqrt((-4-2*so-2*dY*(7+6*dY+so)+u+(2*dY+so+(-1+2*dY)*(3+6*dY+so)*lambda)*u)^2-8*(1+2*dY)*(-2*(1+dY)*(1+2*dY)+(-1+2*dY)*(1+2*dY+so)*lambda*u)*(-1-so+u-lambda*u+2*dY*(-1+lambda*u)))))
      df[length(df$dY)+1,]=c(dY,lambda,so,"Xf",1-((1/(4*(1+2*dY)*(-1-so+u-lambda*u+2*dY*(-1+lambda*u))))*(so*(-2+u-lambda*u+2*dY*(-1+lambda*u))+(1+2*dY)*(-4+u-3*lambda*u+6*dY*(-1+lambda*u))-sqrt((-4-2*so-2*dY*(7+6*dY+so)+u+(2*dY+so+(-1+2*dY)*(3+6*dY+so)*lambda)*u)^2-8*(1+2*dY)*(-2*(1+dY)*(1+2*dY)+(-1+2*dY)*(1+2*dY+so)*lambda*u)*(-1-so+u-lambda*u+2*dY*(-1+lambda*u))))))
    }
  }
}

for(i in c(1,2,3,5)) {df[,i]=as.numeric(df[,i])}

pdf("lambda_OXf_eq.pdf",height=6,width=7)
ggplot(df,aes(x=dY,y=lambda,fill=Eq))+geom_tile()+theme_cowplot()+scale_y_continuous(trans='log10')+
  scale_fill_gradient(low="white",high="red",na.value="gray")+facet_wrap(Chrom~so,ncol=4)
dev.off()

###Figure 3C
df=data.frame(matrix(NA,ncol=4,nrow=0))
names(df)=c("dY","lambda","so","reduction")

u= 0.001; 
for(so in c(0.1,0.2,0.4,0.8)) {
  for(dY in seq(0.05,0.5,0.01)) {
    for(lambda in 10^seq(0,3,0.1)) {
      
      df[length(df$dY)+1,]=c(dY,lambda,so,(1-1/(4*(1+2*dY)*(-1-so+u-lambda*u+2*dY*(-1+lambda*u)))*(so*(-2+u-lambda*u+2*dY*(-1+lambda*u))+(1+2*dY)*(-4+u-3*lambda*u+6*dY*(-1+lambda*u))-sqrt((-4-2*so-2*dY*(7+6*dY+so)+u+(2*dY+so+(-1+2*dY)*(3+6*dY+so)*lambda)*u)^2-8*(1+2*dY)*(-2*(1+dY)*(1+2*dY)+(-1+2*dY)*(1+2*dY+so)*lambda*u)*(-1-so+u-lambda*u+2*dY*(-1+lambda*u)))))/(1-(1+dY)/(1+2*dY)))
      }
  }
}

pdf("lambda_SR_Reduction.pdf",height=4,width=7)
ggplot(df,aes(x=dY,y=lambda,z=1-reduction,fill=1-reduction))+geom_tile()+scale_y_continuous(trans='log10')+theme_cowplot()+
  scale_fill_gradient(low="white",high="red",na.value="gray")+facet_wrap(~so,ncol=4)+geom_contour(col="darkgray",na.rm=TRUE,breaks=c(0.05,seq(0,1,0.2)))+
  geom_text_contour(aes(z=1-reduction),breaks=c(0.05,seq(0,1,0.2)),col="darkgray",nudge_x=0.03,nudge_y=-0.02,rotate=FALSE)+theme(legend.position=c(0.8,0.2))
dev.off()

library(ggplot2)
library(cowplot)
library(reshape2)

##Figure 4 - stability

h=0.0
dO=-0.5

df=data.frame(matrix(NA,ncol=4,nrow=0))
names(df)=c("ssr","t","XfSR","O")

j=1
for(dY in c(0.45)) {
  for(so in c(0.5)) {
    for(ssro in c(0)) {
      for(ssry in c(0.1)) {
        for(ssr in seq(0.1,0.9,0.4)) {
          for(u in c(0.001)) {
            for(gamma in c(10^1)) {
              
              xf=0.999
              xm=0.999
              y=1-u/so
              
              t=0
              
              xflist=c(xf)
              xmlist=c(xm)
              ylist=c(y)
              
              while(t<2000) {
                
                wbarxf=xf*xm+(1-h*ssr)*(xm*(1-xf)+xf*(1-xm))+(1-ssr)*(1-xf)*(1-xm);
                wbarxm=1/2*xf*y+1/2*(1-so)*xf*(1-y)+(1/2+dY)*y*(1-xf)*(1-ssry)+(1/2+dO)*(1-ssro)*(1-y)*(1-xf);
                wbary=1/2*xf*y+1/2*(1-so)*xf*(1-y)+(1/2-dY)*y*(1-xf)*(1-ssry)+(1/2-dO)*(1-ssro)*(1-y)*(1-xf);
                
                xfnext=(xf*xm+1/2*(1-h*ssr)*(xm*(1-xf)+xf*(1-xm)))/wbarxf
                xmnext=(1/2*xf*y+1/2*(1-so)*xf*(1-y))/wbarxm
                ynext=(1/2*xf*y*(1-u)+(1/2-dY)*(1-xf)*y*(1-u*gamma)*(1-ssry))/wbary
                
                xf=xfnext
                xm=xmnext
                y=ynext
                
                xflist=append(xflist,xf)
                xmlist=append(xmlist,xm)
                ylist=append(ylist,y)
                df[length(df$ssr)+1,]=c(ssr,t,1-xf,1-y)
                t=t+1
              }
              
              
              xfsrlist=1-xflist
              olist=1-ylist
              xmsrlist=1-xmlist
              
              j=j+1
              
            }
          }
        }
      }
    }
  }
}

df2=melt(df,id.vars = c("ssr","t"),measure.vars = c("XfSR","O"),variable.name = "Chromosome",value.name = "Equilibrium")

# Figure 4A
pdf("Stability_XfO.pdf")
ggplot(df2[which(df2$ssr %in% c(0.1,0.5,0.9)),],aes(x=t,y=Equilibrium,by=as.factor(ssr)))+geom_line(col=gray(0.5,0.5),aes(linetype=as.factor(ssr)))+
  theme_cowplot()+xlim(0,1000)+theme(legend.position=c(0.2,0.9))+facet_wrap(~Chromosome,nrow=1)+scale_linetype_manual(values=c("dotted","twodash", "solid"))
dev.off()

# Figure 4B
pdf("Stability_XfO2.pdf")
ggplot(df[which(df$ssr %in% c(0.1,0.5,0.9)),],aes(x=O,y=XfSR,by=as.factor(ssr)))+geom_path(col=gray(0.5,0.5))+theme_cowplot()+
  ylim(0,0.6)+xlim(0,1)+
  annotate(geom = 'text', label = 'ssr=0.1', x = 0.75, y = 0.55, hjust = 0, vjust = 1)+
  annotate(geom = 'text', label = 'ssr=0.5', x = 0.25, y = 0.37, hjust = 0, vjust = 1)+
  annotate(geom = 'text', label = 'ssr=0.9', x = 0.03, y = 0.24, hjust = 0, vjust = 1)
dev.off()

## Figure S.1 - stability

d<-read.csv("~/Google Drive/UncklessLab/SRODynamics/Simulations_and_Scripts/gamma_model/nondisjunction.sims.csv",header=FALSE)
names(d)=c("t","u","gamma","ggamm*u","dY","dO","ssr","ssry","ssro","so","h","XfSR","XmSR","O","O-MSB",
           "OStability","XmStability","XfStability","OInvade","XmInvade","XfInvade")

d$XfSR2=d$XfSR
d$XfSR2[which(d$XfStability==0)]=NA
pdf("Stability_Params.pdf")
ggplot(d[intersect(which(d$ssry==0.1),intersect(which(d$dY==0.45),which(d$u==0.001))),],
       aes(x=ssr,y=so,z=XfSR2,fill=XfSR2))+geom_tile()+facet_wrap(gamma~ssro)+theme_cowplot()+
  scale_fill_gradient(low="white",high="red",na.value="gray")+annotate(geom="text",label="Oscillations",srt=90,x=0.2,y=0.8)
dev.off()

### Figure S.2 - 

d<-read.csv("nondisjunction.sims.Yextinction.csv",header=FALSE)

names(d)=c("t","u","gamma","u*gamma","dY","dO","ssr","ssry","ssro","so","h",
           "XfSR","XmSR","O","u/so","Ostable","XmStable","XfSRStable","OInvade","XmInvade",
           "XfInvade")

pdf("SRO_final_frequencies.pdf")
ggplot(d[which(d$Ostable==1),],aes(x=O,y=XfSR))+geom_point(col="darkred")+theme_cowplot()+ylim(0,1)+xlim(0,1)
dev.off()




