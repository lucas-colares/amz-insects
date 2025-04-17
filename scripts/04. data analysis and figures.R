# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # Data analysis script for Colares et al., 2025 # # # # #
# # # # # # # # # # # # # # By L. Colares # # # # # # # # # # # #
source("scripts/00. setup.R")

read.csv("datasets/csv/matrix of insects.csv",h=T,row.names = 1)->insects
read.csv("datasets/csv/landscape metric.csv",h=T,row.names = 1)->land
land = land[match(unique(insects$trap)[order(unique(insects$trap))],land$Trap),]

##### 1. Calculate richness, abundance and composition ----
table(insects$trap,insects$class)->ins_mat
ins_occ=ins_mat
ins_occ[ins_occ>0]<-1
data.frame(class=colnames(ins_mat),hab=c("T","T","A","T","T","A","A"))->habs

vegdist(ins_mat,method = "bray")->d
cmdscale(d,k=2,eig = F)->pcoa
colnames(pcoa)=c("PCoA1","PCoA2")
cmdscale(vegdist(ins_occ,"jaccard"),k=2,eig = F)->pcoa_occ
colnames(pcoa_occ)=c("PCoA1_occ","PCoA2_occ")

full_data = data.frame(Trap=rownames(ins_mat),
                       Long=land$X,
                       Lat=land$Y,
                       hab=c(rep("T",nrow(ins_mat)),rep("A",nrow(ins_mat))),
                       Abun=c(rowSums(ins_mat[,habs$hab=="T"]),rowSums(ins_mat[,habs$hab=="A"])),
                       S=c(rowSums(ins_occ[,habs$hab=="T"]),rowSums(ins_occ[,habs$hab=="A"])),
                       pcoa,
                       pcoa_occ,
                       prop=c(rowSums(ins_mat[,habs$hab=="T"])/rowSums(ins_mat),rowSums(ins_mat[,habs$hab=="A"])/rowSums(ins_mat)),
                       propOcc=c(rowSums(ins_occ[,habs$hab=="T"])/rowSums(ins_occ),rowSums(ins_occ[,habs$hab=="A"])/rowSums(ins_occ)))

full_data$S_n=c(vegan::rarefy(ins_mat[,habs$hab=="T"],sample= min(rowSums(ins_mat))),
                vegan::rarefy(ins_mat[,habs$hab=="A"],sample= min(rowSums(ins_mat))))


full_data = full_data[(full_data$Abun>=4)&(full_data$Abun > full_data$S),]
land_fil = land[match(full_data$Trap,land$Trap),]

full_data$ENS_n = map_dbl(full_data$S_n, ENS_fun, n=4)
full_data$ENS_N = map2_dbl(full_data$S,full_data$Abun, ENS_fun)

full_data$SAD_effect = full_data$ENS_n
full_data$N_effect = full_data$ENS_N-full_data$ENS_n

##### 2. Scale of effect ----
##### 2.1. Richness' SoE ----
c(100,250,500,750,1000,1500,2000,2500,3000)->scales

GAM_S_res=pblapply(scales,function(x){
  gam_data<-data.frame(full_data,Area=land_fil[,paste0("Area",sprintf("%04d",x))])
  
  AICs<-lapply(1:20,function(y){
    gam(SAD_effect~Area*hab+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_S1
    gam(N_effect~Area*hab+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_S2
    return(data.frame(SAD=AIC(GAM_S1),N=AIC(GAM_S2),k=y))
  })
  do.call("rbind",AICs)->k_data
  k_data[order(k_data$SAD),][1,]$k->temp_k
  gam(SAD_effect~Area*hab+s(Long, Lat,k=temp_k),family = "gaussian",data = gam_data)->GAM_S1
  gam(N_effect~Area*hab+s(Long, Lat,k=temp_k),family = "gaussian",data = gam_data)->GAM_S2
  summary(GAM_S1)
  summary(GAM_S2)
  dist(gam_data[,c("Lat","Long")])->distzzz1
  dist(na.omit(gam_data)[,c("Lat","Long")])->distzzz2
  Moran.I(GAM_S1$residuals,as.matrix(distzzz1))->MoranI1
  Moran.I(GAM_S2$residuals,as.matrix(distzzz2))->MoranI2
  
  coefs_mat = rbind(data.frame(Estimate=c(summary(GAM_S1)[[1]],NA),error=c(summary(GAM_S1)[[22]][,2],NA),k=temp_k,stat=c(summary(GAM_S1)[[3]],summary(GAM_S1)[[7]]),P=c(summary(GAM_S1)[[4]],summary(GAM_S1)[[8]]),R2=summary(GAM_S1)[[10]],AIC=AIC(GAM_S1),scale=x,var=c(names(summary(GAM_S1)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_S1)[[1]])),"Random effects"),component="SAD"),
        data.frame(Estimate=c(summary(GAM_S2)[[1]],NA),error=c(summary(GAM_S2)[[22]][,2],NA),k=temp_k,stat=c(summary(GAM_S2)[[3]],summary(GAM_S2)[[7]]),P=c(summary(GAM_S2)[[4]],summary(GAM_S2)[[8]]),R2=summary(GAM_S2)[[10]],AIC=AIC(GAM_S2),scale=x,var=c(names(summary(GAM_S2)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_S2)[[1]])),"Random effects"),component="N"))
  
  random_moran= data.frame(observed=c(MoranI1$observed,MoranI2$observed),
             expected=c(MoranI1$expected,MoranI2$expected),
             sd=c(MoranI1$sd,MoranI2$sd),
             p=c(MoranI1$p.value,MoranI2$p.value),scale=x,
             component=c("SAD","N"))
  
  return(list(moran=random_moran,
              coefs=coefs_mat))
})
do.call("rbind",lapply(GAM_S_res,"[[",1))->moran_S
do.call("rbind",lapply(GAM_S_res,"[[",2))->SoE_S
add_significance(SoE_S,p.col="P",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 0.1,1),
                 symbols = c("****", "***", "**", "*", ".","ns"))->SoE_S

SoE_S[SoE_S$component=="SAD",][order(SoE_S[SoE_S$component=="SAD",]$R2,decreasing = T),][1:5,]
SoE_S[SoE_S$component=="N",][order(SoE_S[SoE_S$component=="N",]$R2,decreasing = T),][1:5,]

SoE_S[SoE_S$component=="SAD",][order(SoE_S[SoE_S$component=="SAD",]$R2,decreasing = T),][1:5,]$scale->sel_SAD
SoE_S[SoE_S$component=="N",][order(SoE_S[SoE_S$component=="N",]$R2,decreasing = T),][1:5,]->sel_N

SoE_S[SoE_S$component=="SAD",]->only_SAD
SoE_S[SoE_S$component=="N",]->only_N

only_SAD[(only_SAD$var=="Area"),]$Estimate/only_SAD[(only_SAD$scale==sel_SAD)&(only_SAD$var=="Area"),]$Estimate
only_N[(only_N$var=="Area"),]$Estimate/only_N[(only_N$scale==sel_SAD)&(only_N$var=="Area"),]$Estimate
only_N[(only_N$var=="Area:habT"),]$Estimate/only_N[(only_N$scale==sel_SAD)&(only_N$var=="Area:habT"),]$Estimate

SoE_S[grepl("Area",rownames(SoE_S)),]->area_coefs_S

rbind(area_coefs_S[area_coefs_S$component=="SAD",][order(area_coefs_S[area_coefs_S$component=="SAD",]$R2,decreasing = T),][1:2,],
area_coefs_S[area_coefs_S$component=="N",][order(area_coefs_S[area_coefs_S$component=="N",]$R2,decreasing = T),][1:2,])->sel_scas
# 
# area_coefs_S$component = gsub("SAD","(c) SAD-component:",area_coefs_S$component)
# area_coefs_S$component = gsub("N","(d) N-component:",area_coefs_S$component)
# sel_scas$component = gsub("SAD","(c) SAD-component:",sel_scas$component)
# sel_scas$component = gsub("N","(d) N-component:",sel_scas$component)

area_coefs_S$hab<-c("A","T")
area_coefs_S$R2[area_coefs_S$component=="SAD"]<-as.numeric(decostand(area_coefs_S$R2[area_coefs_S$component=="SAD"],"range"))
area_coefs_S$R2[area_coefs_S$component=="N"]<-as.numeric(decostand(area_coefs_S$R2[area_coefs_S$component=="N"],"range"))

ScaEffSAD<-ggplot()+
  geom_point(data=area_coefs_S[area_coefs_S$component=="SAD",],aes(y=Estimate,x=scale, size=R2, color=hab))+
  geom_linerange(data=area_coefs_S[area_coefs_S$component=="SAD",],aes(ymin = Estimate-(1.96*error), ymax = Estimate+(1.96*error),x=scale, color=hab))+
  geom_line(data=area_coefs_S[area_coefs_S$component=="SAD",],aes(y=Estimate,x=scale, color=hab),linetype="solid")+
  geom_vline(data=sel_scas[sel_scas$component=="SAD",],aes(xintercept = scale),linetype="dashed",alpha=0.5)+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.25)+
  facet_wrap(.~"(c) SAD-component (SoE):")+
  theme_minimal()+
  scale_color_manual(values=rep(c("#1f96c4","#20684a"),2),labels=c("Aquatic","Terrestrial"))+
  labs(x="Spatial scale (m)",y="Effect size",size="Adj-R²",color="Life cycle")+
  scale_size_continuous(labels=c("Min","Mid-low","Mid","Mid-high","Max"),breaks = seq(min(area_coefs_S$R2),max(area_coefs_S$R2),length.out=5),limits = c(min(area_coefs_S$R2),max(area_coefs_S$R2)))+
  # guides(color="none",size="none")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));ScaEffSAD

ScaEffN<-ggplot()+
  geom_point(data=area_coefs_S[area_coefs_S$component=="N",],aes(y=Estimate,x=scale, size=R2, color=hab))+
  geom_linerange(data=area_coefs_S[area_coefs_S$component=="N",],aes(ymin = Estimate-(1.96*error), ymax = Estimate+(1.96*error),x=scale, color=hab))+
  geom_line(data=area_coefs_S[area_coefs_S$component=="N",],aes(y=Estimate,x=scale, color=hab),linetype="solid")+
  geom_vline(data=sel_scas[sel_scas$component=="N",],aes(xintercept = scale),linetype="dashed",alpha=0.5)+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.25)+
  facet_wrap(.~"(d) N-component (SoE):")+
  theme_minimal()+
  scale_size_continuous(labels=c("Min","Mid-low","Mid","Mid-high","Max"),breaks = seq(min(area_coefs_S$R2),max(area_coefs_S$R2),length.out=5),limits = c(min(area_coefs_S$R2),max(area_coefs_S$R2)))+
  scale_color_manual(values=rep(c("#1f96c4","#20684a"),2),labels=c("Aquatic","Terrestrial"))+
  labs(x="Spatial scale (m)",y="Effect size",size="Adj-R²",color="Life cycle")+
  # guides(color="none",size="none")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));ScaEffN

SADPlt<-ggplot()+
  geom_point(data=full_data,aes(x=land_fil$Area0750,y=log10(SAD_effect),color=hab))+
  geom_smooth(data=full_data,aes(x=land_fil$Area0750,y=log10(SAD_effect),color=hab),method="glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~"(a) SAD-component:")+
  scale_color_manual(values=c("#1f96c4","#20684a"),labels=c("Aquatic","Terrestrial"))+
  theme_minimal()+
  labs(x="Proportion of forest",y="Effective number of groups",color="Life cycle")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));SADPlt

NPlt<-ggplot()+
  geom_point(data=full_data,aes(x=land_fil$Area0100,y=N_effect,color=hab))+
  geom_smooth(data=full_data,aes(x=land_fil$Area0100,y=N_effect,color=hab),method="glm",method.args=c(family="gaussian"))+
  facet_wrap(.~"(b) N-component:")+
  # scale_color_manual(values=c("#0766AD","#6EC207"))+
  scale_color_manual(values=c("#1f96c4","#20684a"))+
  theme_minimal()+
  labs(x="Proportion of forest",y="Effective number of groups",color="Life cycle")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));NPlt

# s_plots<-ggarrange(SPlt,ScaEffS,nrow = 2,heights = c(1,0.6));s_plots

##### 2.2. Composition SoE ----
##### 2.2.1. Coverage-based rarefaction
tagC_T<-min(beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))
tagC_A<-min(beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))

C_stand=round(min(tagC_T, tagC_A)- 0.01, 2)

T_betaC<-beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)
A_betaC<-beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)

box_data<-data.frame(beta=c(T_betaC,A_betaC),group=c(rep("T",9999),rep("A",9999)),type="(a) Coverage-based rarefaction:",var="abun")

T_betaTrue<-betaC:::beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_true"),summarise = F,resamples = 9999)
A_betaTrue<-betaC:::beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_true"),summarise = F,resamples = 9999)

box_data2<-data.frame(beta=c(T_betaTrue,A_betaTrue),group=c(rep("T",9999),rep("A",9999)),type="(b) Whittaker's:",var="abun")

tagC_T<-min(beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))
tagC_A<-min(beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))

C_stand=round(min(tagC_T, tagC_A)- 0.01, 2)

T_betaC<-beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)
A_betaC<-beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)

box_data3<-data.frame(beta=c(T_betaC,A_betaC),group=c(rep("T",9999),rep("A",9999)),type="(a) Coverage-based rarefaction:",var="occ")

T_betaTrue<-betaC:::beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_true"),summarise = F,resamples = 9999)
A_betaTrue<-betaC:::beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_true"),summarise = F,resamples = 9999)

box_data4<-data.frame(beta=c(T_betaTrue,A_betaTrue),group=c(rep("T",9999),rep("A",9999)),type="(b) Whittaker's:",var="occ")

final_box <- rbind(box_data,box_data2,box_data3,box_data4)
gsub("\\(b\\)","(d)",final_box[final_box$var=="abun",]$type)->final_box[final_box$var=="abun",]$type
gsub("\\(a\\)","(c)",final_box[final_box$var=="abun",]$type)->final_box[final_box$var=="abun",]$type

betaC_plt<-ggplot(data = final_box[final_box$var=="abun",],mapping = aes(x=group,y=beta))+
  geom_violin(mapping = aes(color=group,fill=group),trim = F,width=1.5)+
  stat_summary(fun.data=mean_sd,geom="pointrange",color="white")+
  facet_wrap(.~type)+
  scale_color_manual(values=c("#1f96c4","#20684a"))+
  scale_fill_manual(values=c("#1f96c4","#20684a"))+
  scale_x_discrete(labels=c("Aquatic","Terrestrial"))+
  #stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns = T)+ 
  theme_minimal()+
  guides(fill="none",color="none")+
  labs(x="Life cycle",y="β-diversity",color="Life cycle")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));betaC_plt

# Cs1 = lapply(c("F","A|^M","TF"),function(x){
#   Cs2 = lapply(c("T","A"), function(y){
#     sel_hab = ins_mat[rowSums(ins_mat[,habs$hab==y])>0,habs$hab==y]
#     sel_hab = sel_hab[grepl(paste0("^",x),rownames(sel_hab)),]
#     c_targ<-min(beta_stand(sel_hab, setsize = 20, list("C_target"),summarise = F,resamples = 9999))
#     return(c_targ)
#   })
#   return(unlist(Cs2))
# })
# 
# C_stand<-round(min(unlist(Cs1))-0.01,2)
# 
# betaC1 = lapply(c("F","A|^M","TF"),function(x){
#   betaC2 = lapply(c("T","A"), function(y){
#     sel_hab = ins_mat[rowSums(ins_mat[,habs$hab==y])>0,habs$hab==y]
#     sel_hab = sel_hab[grepl(paste0("^",x),rownames(sel_hab)),]
#     T_betaC<-beta_stand(sel_hab, setsize = 20, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)
#     T_betaW<-beta_stand(sel_hab, setsize = 20, list("beta_true"),summarise = F,resamples = 9999)
#     return(data.frame(betaC=T_betaC,betaW=T_betaW,life=y,area=x))
#   })
#   return(do.call("rbind",betaC2))
# })
# do.call("rbind",betaC1)->betaC_est
# betaC_est<-melt(betaC_est)
# 
# betaC_plt2<-ggplot(data = betaC_est,mapping = aes(x=area,y=value))+
#   geom_violin(mapping = aes(color=life,fill=life),trim = F,width=1.5)+
#   stat_summary(fun.data=mean_ci, mult=1,geom="pointrange",color="white")+
#   facet_wrap(.~variable)+
#   scale_color_manual(values=c("#1f96c4","#20684a"))+
#   scale_fill_manual(values=c("#1f96c4","#20684a"))+
#   #scale_x_discrete(labels=c("Aquatic","Terrestrial"))+
#   #stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns = T)+ 
#   theme_minimal()+
#   guides(fill="none",color="none")+
#   labs(x="Life cycle",y="β-diversity",color="Life cycle")+
#   theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));betaC_plt2

##### 2.2.2. PCoA Axis
cmdscale(d,k=2,eig = T,add=T)->pcoa2

cmdscale(vegdist(ins_occ,"jaccard"),k=2,eig = T,add=T)->pcoa_occ

GAM_Comp_res=pblapply(scales,function(x){
  gam_data<-data.frame(full_data[full_data$hab=="T",],Area=land_fil[full_data$hab=="T",paste0("Area",sprintf("%04d",x))])
  
  AICs<-lapply(1:20,function(y){
    gam(PCoA1~Area+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_Comp
    gam(PCoA1_occ~Area+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_Comp2
    return(c(AIC(GAM_Comp,GAM_Comp2)))
  })
  
  data.frame(k=1:20,AIC1=unlist(lapply(AICs,"[[",1)),AIC2=unlist(lapply(AICs,"[[",2)))->k_data
  k_data[order(k_data$AIC1),][1,]$k->temp_k1
  k_data[order(k_data$AIC2),][1,]$k->temp_k2
  gam(PCoA1~Area+s(Long, Lat,k=temp_k1),family = "gaussian",data = gam_data)->GAM_Comp
  gam(PCoA1_occ~Area+s(Long, Lat,k=temp_k2),family = "gaussian",data = gam_data)->GAM_Comp2
  dist(gam_data[,c("Lat","Long")])->distzzz
  Moran.I(GAM_Comp$residuals,as.matrix(distzzz))->MoranI1
  Moran.I(GAM_Comp2$residuals,as.matrix(distzzz))->MoranI2
  
  full_coefs = data.frame(Estimate=c(summary(GAM_Comp)[[1]],NA,summary(GAM_Comp2)[[1]],NA),
             k=c(rep(temp_k1,((length(summary(GAM_Comp)[[1]])+1))),rep(temp_k2,((length(summary(GAM_Comp)[[1]])+1)))),
             stat=c(summary(GAM_Comp)[[3]],summary(GAM_Comp)[[7]],summary(GAM_Comp2)[[3]],summary(GAM_Comp2)[[7]]),
             P=c(summary(GAM_Comp)[[4]],summary(GAM_Comp)[[8]],summary(GAM_Comp2)[[4]],summary(GAM_Comp2)[[8]]),
             R2=c(rep(summary(GAM_Comp)[[10]],((length(summary(GAM_Comp)[[1]])+1))),rep(summary(GAM_Comp2)[[10]],((length(summary(GAM_Comp)[[1]])+1)))),
             AIC=c(rep(AIC(GAM_Comp),((length(summary(GAM_Comp)[[1]])+1))),rep(AIC(GAM_Comp2),((length(summary(GAM_Comp)[[1]])+1)))),
             scale=x,
             var=c(names(summary(GAM_Comp)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_Comp)[[1]])),"Random effects"),
             data=c(rep("Abundance",((length(summary(GAM_Comp)[[1]])+1))),rep("Occurrence",((length(summary(GAM_Comp)[[1]])+1)))))
  
  return(list(moran=data.frame(observed=c(MoranI1$observed,MoranI2$observed),expected=c(MoranI1$expected,MoranI2$expected),sd=c(MoranI1$sd,MoranI2$sd),p=c(MoranI1$p.value,MoranI2$p.value),scale=x,data=c("Abundance","Occurrence")),
              coefs=full_coefs))
})
do.call("rbind",lapply(GAM_Comp_res,"[[",1))->moran_Comp
do.call("rbind",lapply(GAM_Comp_res,"[[",2))->SoE_Comp
add_significance(SoE_Comp,p.col="P",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 0.1,1),
                 symbols = c("****", "***", "**", "*", ".","ns"))->SoE_Comp

SoE_Comp[SoE_Comp$data=="Abundance",][order(SoE_Comp[SoE_Comp$data=="Abundance",]$R2,decreasing = T),][1:3,]
SoE_Comp[SoE_Comp$data=="Occurrence",][order(SoE_Comp[SoE_Comp$data=="Occurrence",]$R2,decreasing = T),][1:3,]

unique(SoE_Comp[SoE_Comp$data=="Abundance",][order(SoE_Comp[SoE_Comp$data=="Abundance",]$R2,decreasing = T),][1:3,]$scale)->sel_Abun
unique(SoE_Comp[SoE_Comp$data=="Occurrence",][order(SoE_Comp[SoE_Comp$data=="Occurrence",]$R2,decreasing = T),][1:3,]$scale)->sel_Occ

SoE_Comp[SoE_Comp$data=="Abundance",]->only_Abun
SoE_Comp[SoE_Comp$data=="Occurrence",]->only_Occ

only_Abun[(only_Abun$var=="Area"),]$Estimate/only_Abun[(only_Abun$scale==sel_Abun)&(only_Abun$var=="Area"),]$Estimate
only_Occ[(only_Occ$var=="Area"),]$Estimate/only_Occ[(only_Occ$scale==sel_Occ)&(only_Occ$var=="Area"),]$Estimate

SoE_Comp[grepl("Area",SoE_Comp$var),]->area_coefs_comp

# area_coefs_Comp$hab<-c("(e) Compdance's SoE","(f) Compdance's SoE")

for(x in unique(area_coefs_comp$data)){
  area_coefs_comp[area_coefs_comp$data==x,]$R2<-vegan::decostand(area_coefs_comp[area_coefs_comp$data==x,]$R2,method = "range")
}

ScaEffComp<-ggplot()+
  geom_point(data=area_coefs_comp,aes(y=Estimate,x=scale, size=R2,color=data))+
  geom_line(data=area_coefs_comp,aes(y=Estimate,x=scale,color=data),linetype="solid")+
  geom_vline(data=area_coefs_comp[order(area_coefs_comp$R2,decreasing = T),][1,],aes(xintercept = scale),linetype="dashed")+
  theme_minimal()+
  facet_wrap(.~"(e) Composition's SoE (PCoA):")+
  MoMAColors::scale_color_moma_d(palette_name = "VanGogh")+
  scale_y_continuous(breaks=seq(0,5,1.25))+
  scale_size_continuous(labels=c("Min","Mid-low","Mid","Mid-high","Max"),breaks = seq(min(area_coefs_comp$R2),max(area_coefs_comp$R2),length.out=5),limits = c(min(area_coefs_comp$R2),max(area_coefs_comp$R2)))+
  labs(x="Spatial scale (meters)",y="Effect size",size="Adj-R²",color="Biological\ndata")+
  theme(legend.justification = "left",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=9), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));ScaEffComp

CompPlt1<-ggplot()+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.4)+
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.4)+
  #geom_point(data=comp_data,aes(x=PCoA1,y=PCoA2,fill=Area),shape=21,color="black")+
  geom_point(data=full_data[full_data$hab=="T",],aes(x=PCoA1,y=PCoA2,color=prop,size=land_fil$Area0100[full_data$hab=="T"]))+
  theme_minimal()+
  facet_wrap(.~"(b) Abundance change (PCoA):")+
  # guides(color = guide_colorbar(barwidth = 5.5,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"),
  #        size=guide_legend(barwidth = 1,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"))+
  # #guides(fill=guide_legend(ncol=2),size=guide_legend(ncol=2))+
  #scale_color_gradient(low = "#0766AD",high = "#6EC207",guide = "colorbar")+
  scale_color_gradient(low = "#1f96c4",high = "#20684a",guide = "colorbar")+
  labs(x=paste0("PCoA 1 (",round((pcoa2$eig[1]/sum(pcoa2$eig))*100,1),"%)"),y=paste0("PCoA 2 (",round((pcoa2$eig[2]/sum(pcoa2$eig))*100,1),"%)"),color="Proportion\nof terrestrial\nindividuals",size="Proportion\nof forest")+
  theme(legend.justification = "left",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=9), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));CompPlt1

CompPlt2<-ggplot()+
  #geom_point(data=comp_data,aes(x=PCoA1,y=PCoA2,fill=Area),shape=21,color="black")+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.4)+
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.4)+
  geom_point(data=full_data[full_data$hab=="T",],aes(x=PCoA1_occ,y=PCoA2_occ,color=propOcc,size=land_fil$Area0100[full_data$hab=="T"]))+
  theme_minimal()+
  facet_wrap(.~"(a) Compositional change (PCoA):")+
  # guides(color = guide_colorbar(barwidth = 5.5,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"),
  #        size=guide_legend(barwidth = 1,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"))+
  # #guides(fill=guide_legend(ncol=2),size=guide_legend(ncol=2))+
  #scale_color_gradient(low = "#0766AD",high = "#6EC207",guide = "colorbar")+
  scale_color_gradient(low = "#1f96c4",high = "#20684a",guide = "colorbar")+
  labs(x=paste0("PCoA 1 (",round((pcoa3$eig[1]/sum(pcoa3$eig))*100,1),"%)"),y=paste0("PCoA 2 (",round((pcoa3$eig[2]/sum(pcoa3$eig))*100,1),"%)"),color="Proportion\nof terrestrial\ngroups",size="Proportion\nof forest")+
  theme(legend.justification = "left",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=9), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));CompPlt2

##### Taxon-specific SoE (Abundance) ----
taxa_mat = melt(ins_mat)
colnames(taxa_mat) = c("Trap","Taxa","Abun")

GAM_taxa_res=pblapply(scales,function(x){
  gam_data<-data.frame(taxa_mat,Long=land$X,Lat=land$Y,Area=land[,paste0("Area",sprintf("%04d",x))])
  
  taxa_gam_res<-lapply(unique(gam_data$Taxa),function(w){
    taxa_gam = gam_data[gam_data$Taxa==w,]
    
    AICs<-lapply(1:20,function(y){
      gam(Abun~Area+s(Long, Lat,k=y),family = "poisson",data = taxa_gam)->GAM_Abun
      return(AIC(GAM_Abun))
    })
    data.frame(k=1:20,AIC=unlist(AICs))->k_data
    k_data[order(k_data$AIC),][1,]$k->temp_k
    gam(Abun~Area+s(Long, Lat,k=temp_k),family = "poisson",data = taxa_gam)->GAM_Abun
    summary(GAM_Abun)
    dist(taxa_gam[,c("Lat","Long")])->distzzz
    Moran.I(GAM_Abun$residuals,as.matrix(distzzz))->MoranI
    
    return(list(moran=data.frame(observed=MoranI$observed,expected=MoranI$expected,sd=MoranI$sd,p=MoranI$p.value,scale=x, taxa=w),
                coefs=data.frame(Estimate=c(summary(GAM_Abun)[[1]],NA),k=temp_k,stat=c(summary(GAM_Abun)[[3]],summary(GAM_Abun)[[7]]),P=c(summary(GAM_Abun)[[4]],summary(GAM_Abun)[[8]]),R2=summary(GAM_Abun)[[10]],AIC=AIC(GAM_Abun),scale=x,var=c(names(summary(GAM_Abun)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_Abun)[[1]])),"Random effects"), taxa=w)))
  })
  do.call("rbind",lapply(taxa_gam_res,"[[",1))->moran_taxa
  do.call("rbind",lapply(taxa_gam_res,"[[",2))->SoE_taxa
  return(list(moran_taxa,SoE_taxa))
})
do.call("rbind",lapply(GAM_taxa_res,"[[",1))->moran_taxa
do.call("rbind",lapply(GAM_taxa_res,"[[",2))->SoE_taxa
add_significance(SoE_taxa,p.col="P",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 0.1,1),
                 symbols = c("****", "***", "**", "*", ".","ns"))->SoE_taxa

area_coefs_abun_ord = SoE_taxa[grepl("Area",rownames(SoE_taxa)),]

area_abun<-{}
for(x in 1:length(unique(area_coefs_abun_ord$taxa))){
  area_coefs_abun_ord[area_coefs_abun_ord$taxa==unique(area_coefs_abun_ord$taxa)[x],]->sel_tax
  area_abun[[length(area_abun)+1]]<-sel_tax[order(sel_tax$R2,decreasing = T),][1,]
}
do.call("rbind",area_abun)->area_abun

eff_abun_soe<-lapply(unique(area_coefs_abun_ord$taxa),function(x){
  area_abun[area_abun$taxa==x,]$scale->sel_scale
  return(area_coefs_abun_ord[area_coefs_abun_ord$taxa==x,]$Estimate/area_coefs_abun_ord[(area_coefs_abun_ord$taxa==x)&(area_coefs_abun_ord$scale==sel_scale),]$Estimate)
})
min(unlist(eff_abun_soe))
(max(unlist(eff_abun_soe))-1)*100

taxa_mat$Area<-NA
for(x in unique(taxa_mat$Taxa)){
  taxa_mat[taxa_mat$Taxa==x,]$Area<-land[,paste0("Area",sprintf("%04d",area_abun[area_abun$taxa==x,]$scale))]
  area_coefs_abun_ord[area_coefs_abun_ord$taxa==x,]$R2<-vegan::decostand(area_coefs_abun_ord[area_coefs_abun_ord$taxa==x,]$R2,method = "range")
}

c("Bees and wasps","Beetles","Cicadas","Flies","Caddisflies","Mayflies","Mosquitoes")->fantasy_names

c("#5ebe96","#36ae7c","#2b8b63","#20684a","#1f96c4","#18789c","#1b63b2")->my_cols

all_names = data.frame(Taxa=unique(taxa_mat$Taxa)[order(unique(taxa_mat$Taxa))],Fantasy=c("Flies","Beetles","Mayflies","Cicadas","Bees and wasps","Mosquitoes","Caddisflies"))
all_names = all_names[order(all_names$Fantasy),]
all_names = all_names[c(1:2,4:5,3,6:7),]
all_names$Format = paste0("(",letters[1:7],") ",all_names$Fantasy,":")

all_names[match(taxa_mat$Taxa,all_names$Taxa),3] -> taxa_mat$Taxa

as.character(area_coefs_abun_ord$taxa)->area_coefs_abun_ord$taxa
area_coefs_abun_ord$taxa = factor(area_coefs_abun_ord$taxa,levels=all_names$Taxa)

as.character(area_abun$taxa)->area_abun$taxa
area_abun$taxa = factor(area_abun$taxa,levels=all_names$Taxa)

ScaEffOrd_abun<-ggplot()+
  geom_point(data=area_coefs_abun_ord,aes(y=Estimate,x=scale, size=R2, fill=taxa,shape=taxa,color=taxa))+
  geom_line(data=area_coefs_abun_ord,aes(y=Estimate,x=scale, color=taxa, group=taxa),linetype="solid")+
  geom_vline(data=area_abun[grepl("Area",rownames(area_abun)),],aes(xintercept = scale,color=taxa),linetype="dashed")+
  facet_wrap(.~"(h) Abundance's SOE:")+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.4)+
  theme_minimal()+
  scale_size_continuous(labels=c("Min","Mid-low","Mid","Mid-high","Max"),breaks = seq(min(area_coefs_abun_ord$R2),max(area_coefs_abun_ord$R2),length.out=5),limits = c(min(area_coefs_abun_ord$R2),max(area_coefs_abun_ord$R2)))+
  guides(size=guide_legend(order = 2),color=guide_legend(order = 1),fill=guide_legend(order = 1),shape=guide_legend(order = 1))+
  scale_color_manual(values=my_cols,labels=fantasy_names)+
  scale_fill_manual(values=my_cols,labels=fantasy_names)+
  scale_shape_manual(values=c(21,22,23,24,25,8,13),labels=fantasy_names)+
  labs(x="Spatial scale (m)",y="Effect size",size="Adj-R²",fill="Taxa",color="Taxa",shape="Taxa")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));ScaEffOrd_abun

plt_data = taxa_mat[taxa_mat$Abun<500,]
plt_data1 = plt_data[grepl("Caddisflies|Mayflies|Mosquitoes",plt_data$Taxa),]
plt_data2 = plt_data[!grepl("Caddisflies|Mayflies|Mosquitoes",plt_data$Taxa),]

AbunOrdPlt1<-ggplot()+
  geom_point(data=plt_data2,aes(x=Area,y=log10(Abun+1),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=plt_data2,aes(x=Area,y=log10(Abun+1),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~Taxa,scales = "free_y",nrow = 1)+
  scale_color_manual(values=my_cols,labels=fantasy_names)+
  scale_fill_manual(values=my_cols,labels=fantasy_names)+
  scale_linetype(labels=fantasy_names)+
  theme_minimal()+
  scale_y_continuous(labels = c(0,10,100),breaks = c(0,1,2))+
  #guides(color="none",fill="none",color="none",shape="none")+
  scale_shape_manual(values=c(21,22,23,24,25,8,13),labels=fantasy_names)+
  labs(x="Proportion of forest",y="Number of individuals",color="Biological\ngroup",fill="Biological\ngroup",shape="Biological\ngroup",linetype="Biological\ngroup")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));AbunOrdPlt1

AbunOrdPlt2<-ggplot()+
  geom_point(data=plt_data1,aes(x=Area,y=log10(Abun+1),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=plt_data1,aes(x=Area,y=log10(Abun+1),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~Taxa,scales = "free_y",nrow = 1)+
  scale_color_manual(values=my_cols[5:7],labels=fantasy_names[5:7])+
  scale_fill_manual(values=my_cols[5:7],labels=fantasy_names[5:7])+
  scale_linetype(labels=fantasy_names[5:7])+
  theme_minimal()+
  scale_y_continuous(labels = c(0,10,100),breaks = c(0,1,2))+
  #guides(color="none",fill="none",color="none",shape="none")+
  scale_shape_manual(values=c(25,8,13),labels=fantasy_names)+
  labs(x="Proportion of forest",y="Number of individuals",color="Biological\ngroup",fill="Biological\ngroup",shape="Biological\ngroup",linetype="Biological\ngroup")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));AbunOrdPlt2

#### Taxon-specific SoE (Body size) ----
# Phylogenetic tree according to http://dx.doi.org/10.1126/science.1257570
(insects$xmax-insects$xmin)*(insects$ymax-insects$ymin)->insects$body_area

cat("(((((((Brachycera:0.2,Nematocera:0.2):0.2,", "Trichoptera:0.4):0.2,Coleoptera:0.6):0.4,Hymenoptera:1.0):0.6,Hemiptera:1.6):0.6,Ephemeroptera:2.2):0.6);", file = "ex.tre", sep = "\n")
tree.ins <- read.tree("ex.tre")
plot(tree.ins)

size_ins = aggregate(insects$body_area,list(insects$class),max)[,2]
names(size_ins) = unique(insects$class)[order(unique(insects$class))]

K_result <- phylosig(tree.ins, size_ins, method = "lambda", test = TRUE,niter = 9999)
print(K_result)

aggregate(insects$body_area,list(Trap=insects$trap,Taxa=insects$class),max)->size_mat
colnames(size_mat) = c("Trap","Taxa","Size")
size_mat[order(size_mat$Taxa,size_mat$Trap),]->size_mat
expand.grid(unique(insects$trap),unique(insects$class))->size_empty
size_empty$size<-NA
size_empty$size<-size_mat[match(paste0(size_empty$Var1,size_empty$Var2),paste0(size_mat$Trap,size_mat$Taxa)),3]
size_mat = size_empty
colnames(size_mat)[1:3] = c("Trap","Taxa","Size")

GAM_taxa_size=pblapply(scales,function(x){
  gam_data<-data.frame(size_mat,Long=land$X,Lat=land$Y,Area=land[,paste0("Area",sprintf("%04d",x))])
  
  taxa_gam_res<-lapply(unique(gam_data$Taxa),function(w){
    taxa_gam = gam_data[gam_data$Taxa==w,]
    taxa_gam = na.omit(taxa_gam)
    
    AICs<-lapply(1:20,function(y){
      gam(Size~Area+s(Long, Lat,k=y),family = "quasipoisson",data = taxa_gam)->GAM_Size
      return(AIC(GAM_Size))
    })
    data.frame(k=1:20,AIC=unlist(AICs))->k_data
    k_data[order(k_data$AIC),][1,]$k->temp_k
    gam(Size~Area+s(Long, Lat,k=temp_k),family = "quasipoisson",data = taxa_gam)->GAM_Size
    summary(GAM_Size)
    dist(taxa_gam[,c("Lat","Long")])->distzzz
    Moran.I(GAM_Size$residuals,as.matrix(distzzz))->MoranI
    
    return(list(moran=data.frame(observed=MoranI$observed,expected=MoranI$expected,sd=MoranI$sd,p=MoranI$p.value,scale=x, taxa=w),
                coefs=data.frame(Estimate=c(summary(GAM_Size)[[1]],NA),k=temp_k,stat=c(summary(GAM_Size)[[3]],summary(GAM_Size)[[7]]),P=c(summary(GAM_Size)[[4]],summary(GAM_Size)[[8]]),R2=summary(GAM_Size)[[10]],AIC=AIC(GAM_Size),scale=x,var=c(names(summary(GAM_Size)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_Size)[[1]])),"Random effects"), taxa=w)))
  })
  do.call("rbind",lapply(taxa_gam_res,"[[",1))->moran_taxa
  do.call("rbind",lapply(taxa_gam_res,"[[",2))->SoE_taxa
  return(list(moran_taxa,SoE_taxa))
})
do.call("rbind",lapply(GAM_taxa_size,"[[",1))->moran_size
do.call("rbind",lapply(GAM_taxa_size,"[[",2))->SoE_size
add_significance(SoE_size,p.col="P",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 0.1,1),
                 symbols = c("****", "***", "**", "*", ".","ns"))->SoE_size

area_coefs_size_ord = SoE_size[grepl("Area",rownames(SoE_size)),]

area_size<-{}
for(x in 1:length(unique(area_coefs_size_ord$taxa))){
  area_coefs_size_ord[area_coefs_size_ord$taxa==unique(area_coefs_size_ord$taxa)[x],]->sel_tax
  area_size[[length(area_size)+1]]<-sel_tax[order(sel_tax$R2,decreasing = T),][1,]
}
do.call("rbind",area_size)->area_size
area_size

size_mat$Area<-NA
for(x in unique(size_mat$Taxa)){
  size_mat[size_mat$Taxa==x,]$Area<-land[,paste0("Area",sprintf("%04d",area_size[area_size$taxa==x,]$scale))]
  area_coefs_size_ord[area_coefs_size_ord$taxa==x,]$R2<-vegan::decostand(area_coefs_size_ord[area_coefs_size_ord$taxa==x,]$R2,method = "range")
}

as.character(size_mat$Taxa)->size_mat$Taxa
size_mat$Taxa = factor(size_mat$Taxa,levels=all_names$Taxa)

as.character(area_coefs_size_ord$taxa)->area_coefs_size_ord$taxa
area_coefs_size_ord$taxa = factor(area_coefs_size_ord$taxa,levels=all_names$Taxa)

as.character(area_size$taxa)->area_size$taxa
area_size$taxa = factor(area_size$taxa,levels=all_names$Taxa)

eff_size_soe<-lapply(unique(area_coefs_size_ord$taxa),function(x){
  area_size[area_size$taxa==x,]$scale->sel_scale
  return(area_coefs_size_ord[area_coefs_size_ord$taxa==x,]$Estimate/area_coefs_size_ord[(area_coefs_size_ord$taxa==x)&(area_coefs_size_ord$scale==sel_scale),]$Estimate)
})
min(unlist(eff_size_soe))
(max(unlist(eff_size_soe))-1)*100

ScaEffOrd_size<-ggplot()+
  geom_point(data=area_coefs_size_ord,aes(y=Estimate,x=scale, size=R2, fill=taxa,shape=taxa,color=taxa))+
  geom_line(data=area_coefs_size_ord,aes(y=Estimate,x=scale, color=taxa, group=taxa),linetype="solid")+
  geom_vline(data=area_size[grepl("Area",rownames(area_size)),],aes(xintercept = scale,color=taxa),linetype="dashed")+
  facet_wrap(.~"(f) Body size's SOE:")+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.4)+
  theme_minimal()+
  scale_size_continuous(labels=c("Min","Mid-low","Mid","Mid-high","Max"),breaks = seq(min(area_coefs_size_ord$R2),max(area_coefs_size_ord$R2),length.out=5),limits = c(min(area_coefs_size_ord$R2),max(area_coefs_size_ord$R2)))+
  guides(size=guide_legend(order = 2),color=guide_legend(order = 1),fill=guide_legend(order = 1),shape=guide_legend(order = 1))+
  scale_color_manual(values=my_cols,labels=fantasy_names)+
  scale_fill_manual(values=my_cols,labels=fantasy_names)+
  scale_shape_manual(values=c(21,22,23,24,25,8,13),labels=fantasy_names)+
  labs(x="Spatial scale (m)",y="Effect size",size="Adj-R²",fill="Taxa",color="Taxa",shape="Taxa")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));ScaEffOrd_size

size_mat = size_mat[!is.na(size_mat$Size),]
size_mat[(size_mat$Taxa!="Trichoptera")&(size_mat$Taxa!="Hemiptera"),]->size_plt
size_plt$Taxa<-all_names[match(size_plt$Taxa,all_names$Taxa),2]
size_plt$Taxa = factor(size_plt$Taxa,levels = c("Bees and wasps","Beetles","Flies","Mayflies","Mosquitoes"))

size_plt$Taxa = paste0("(",letters[1:5],"): ",levels(size_plt$Taxa))[match(size_plt$Taxa,levels(size_plt$Taxa))]

SizeOrdPlt1<-ggplot()+
  geom_point(data=size_plt[!grepl("Mayflies|Mosquitoes",size_plt$Taxa),],aes(x=Area,y=log10(Size),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=size_plt[!grepl("Mayflies|Mosquitoes",size_plt$Taxa),],aes(x=Area,y=log10(Size),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~Taxa,scales = "free_y")+
  scale_color_manual(values=my_cols[-c(4,7)],labels=fantasy_names[-c(4,7)])+
  scale_fill_manual(values=my_cols[-c(4,7)],labels=fantasy_names[-c(4,7)])+
  scale_linetype(labels=fantasy_names[-c(4,7)])+
  theme_minimal()+
  #scale_y_continuous(labels = c(0,10,100),breaks = c(0,1,2))+
  #guides(color="none",fill="none",color="none",shape="none")+
  scale_shape_manual(values=c(21,22,23,24,25,8,13),labels=fantasy_names[-c(4,7)])+
  labs(x="Proportion of forest",y="Area of the bounding box",color="Biological\ngroup",fill="Biological\ngroup",shape="Biological\ngroup",linetype="Biological\ngroup")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));SizeOrdPlt1

SizeOrdPlt2<-ggplot()+
  geom_point(data=size_plt[grepl("Mayflies|Mosquitoes",size_plt$Taxa),],aes(x=Area,y=log10(Size),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=size_plt[grepl("Mayflies|Mosquitoes",size_plt$Taxa),],aes(x=Area,y=log10(Size),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~Taxa,scales = "free_y")+
  scale_color_manual(values=my_cols[-c(1:3,4,7)],labels=fantasy_names[-c(1:3,4,7)])+
  scale_fill_manual(values=my_cols[-c(1:3,4,7)],labels=fantasy_names[-c(1:3,4,7)])+
  scale_linetype(labels=fantasy_names[-c(1:3,4,7)])+
  theme_minimal()+
  #scale_y_continuous(labels = c(0,10,100),breaks = c(0,1,2))+
  #guides(color="none",fill="none",color="none",shape="none")+
  scale_shape_manual(values=c(8,13),labels=fantasy_names[-c(1:3,4,7)])+
  labs(x="Proportion of forest",y="Area of the bounding box",color="Biological\ngroup",fill="Biological\ngroup",shape="Biological\ngroup",linetype="Biological\ngroup")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));SizeOrdPlt2

##### Arrange and save plots ----
arr1 = ggarrange(AbunOrdPlt1,ggarrange(AbunOrdPlt2,ScaEffOrd_abun+theme(legend.position = "none"),widths = c(1,0.4)),nrow = 2,common.legend = T,legend.grob = get_legend(ScaEffOrd_abun),legend = "right"); arr1
ggsave(filename = "figures/abundance_order.tif",plot = arr1,units = "in",width = 13.5,height = 6,dpi = 600)

arr2 = ggarrange(SizeOrdPlt1,ggarrange(SizeOrdPlt2,ScaEffOrd_size+theme(legend.position = "none"),widths = c(1,0.6)),nrow = 2,common.legend = T,legend.grob = get_legend(ScaEffOrd_size),legend = "right"); arr2
ggsave(filename = "figures/size_order.tif",plot = arr2,units = "in",width = 11,height = 6,dpi = 600)

arr3 = ggarrange(SADPlt,NPlt,common.legend = T,legend = "right");arr3
arr4 = ggarrange(ScaEffSAD,ScaEffN,common.legend = T,legend = "right");arr4
arr5 = ggarrange(arr2,arr3,nrow = 2); arr5
ggsave(filename = "figures/SAD+N.tif",plot = arr5,units = "in",width = 8.5,height = 6,dpi = 600)

arr6 = ggarrange(CompPlt2,CompPlt1,betaC_plt,ScaEffComp,nrow = 2,ncol = 2,heights = c(1,0.85),align = "h");arr6
ggsave(filename = "figures/composition_plot.tif",plot = arr6,units = "in",width = 11.75,height = 7.5,dpi = 600)