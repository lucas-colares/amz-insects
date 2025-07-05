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

##### SAD & N ----
gam_data<-data.frame(full_data,Area=land_fil[,"Area0100"])

AICs<-lapply(1:20,function(y){
  gam(SAD_effect~Area+hab+Area:hab+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_S1
  gam(N_effect~Area+hab+Area:hab+s(Long, Lat,k=y),family = "gaussian",data = gam_data)->GAM_S2
  return(data.frame(SAD=AIC(GAM_S1),N=AIC(GAM_S2),k=y))
})
do.call("rbind",AICs)->k_data
k_data[order(k_data$SAD),][1,]$k->temp_k
k_data[order(k_data$N),][1,]$k->temp_k
gam(SAD_effect~Area+hab+Area:hab+s(Long, Lat,k=temp_k),family = "gaussian",data = gam_data)->GAM_S1
gam(N_effect~Area+hab+Area:hab+s(Long, Lat,k=temp_k),family = "gaussian",data = gam_data)->GAM_S2
summary(GAM_S1)
summary(GAM_S2)
dist(gam_data[,c("Lat","Long")])->distzzz1
dist(na.omit(gam_data)[,c("Lat","Long")])->distzzz2
Moran.I(GAM_S1$residuals,as.matrix(distzzz1))->MoranI1
Moran.I(GAM_S2$residuals,as.matrix(distzzz2))->MoranI2

coefs_mat = rbind(data.frame(Estimate=c(summary(GAM_S1)[[1]],NA),error=c(summary(GAM_S1)[[22]][,2],NA),k=temp_k,stat=c(summary(GAM_S1)[[3]],summary(GAM_S1)[[7]]),P=c(summary(GAM_S1)[[4]],summary(GAM_S1)[[8]]),R2=summary(GAM_S1)[[10]],AIC=AIC(GAM_S1),var=c(names(summary(GAM_S1)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_S1)[[1]])),"Random effects"),component="SAD"),
                  data.frame(Estimate=c(summary(GAM_S2)[[1]],NA),error=c(summary(GAM_S2)[[22]][,2],NA),k=temp_k,stat=c(summary(GAM_S2)[[3]],summary(GAM_S2)[[7]]),P=c(summary(GAM_S2)[[4]],summary(GAM_S2)[[8]]),R2=summary(GAM_S2)[[10]],AIC=AIC(GAM_S2),var=c(names(summary(GAM_S2)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_S2)[[1]])),"Random effects"),component="N"))

random_moran= data.frame(observed=c(MoranI1$observed,MoranI2$observed),
                         expected=c(MoranI1$expected,MoranI2$expected),
                         sd=c(MoranI1$sd,MoranI2$sd),
                         p=c(MoranI1$p.value,MoranI2$p.value),
                         component=c("SAD","N"))

coefs_mat
random_moran

SADPlt<-ggplot()+
  geom_point(data=full_data,aes(x=land_fil$Area0100,y=log(SAD_effect),color=hab))+
  geom_smooth(data=full_data,aes(x=land_fil$Area0100,y=log(SAD_effect),color=hab),method="glm",method.args=list(family="quasipoisson"))+
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

##### 2. Composition ----
##### 2.2.1. Coverage-based rarefaction
tagC_T<-min(beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))
tagC_A<-min(beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))

C_stand=round(min(tagC_T, tagC_A)- 0.01, 2)

T_betaC<-beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)
A_betaC<-beta_stand(ins_mat[rowSums(ins_mat[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)

box_data<-data.frame(beta=c(T_betaC,A_betaC),group=c(rep("T",9999),rep("A",9999)),type="(a) Coverage-based rarefaction:",var="abun")

tagC_T<-min(beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))
tagC_A<-min(beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("C_target"),summarise = F,resamples = 9999))

C_stand=round(min(tagC_T, tagC_A)- 0.01, 2)

T_betaC<-beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="T"])>0,habs$hab=="T"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)
A_betaC<-beta_stand(ins_occ[rowSums(ins_occ[,habs$hab=="A"])>0,habs$hab=="A"], setsize = 99, list("beta_C"),args = list(C= C_stand),summarise = F,resamples = 9999)

box_data3<-data.frame(beta=c(T_betaC,A_betaC),group=c(rep("T",9999),rep("A",9999)),type="(a) Coverage-based rarefaction:",var="occ")

final_box <- rbind(box_data,box_data3)

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

##### 2.2. PCoA Axis
cmdscale(d,k=2,eig = T,add=T)->pcoa2

cmdscale(vegdist(ins_occ,"jaccard"),k=2,eig = T,add=T)->pcoa_occ

gam_data<-data.frame(full_data[full_data$hab=="T",],Area=land_fil[full_data$hab=="T","Area0100"])

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

GAM_Comp_res<-(list(moran=data.frame(observed=c(MoranI1$observed,MoranI2$observed),expected=c(MoranI1$expected,MoranI2$expected),sd=c(MoranI1$sd,MoranI2$sd),p=c(MoranI1$p.value,MoranI2$p.value),scale=x,data=c("Abundance","Occurrence")),
            coefs=full_coefs))

GAM_Comp_res$moran->moran_Comp
GAM_Comp_res$coefs->SoE_Comp

CompPlt1<-ggplot()+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.4)+
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.4)+
  #geom_point(data=comp_data,aes(x=PCoA1,y=PCoA2,fill=Area),shape=21,color="black")+
  geom_point(data=full_data[full_data$hab=="T",],aes(x=PCoA1,y=PCoA2,color=prop,size=land_fil$Area0100[full_data$hab=="T"]))+
  theme_minimal()+
  facet_wrap(.~"(c) Abundance change (PCoA):")+
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
  facet_wrap(.~"(b) Compositional change (PCoA):")+
  # guides(color = guide_colorbar(barwidth = 5.5,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"),
  #        size=guide_legend(barwidth = 1,barheight = 1,label.position = "bottom",title.position = "top", title.vjust = 1,direction = "vertical"))+
  # #guides(fill=guide_legend(ncol=2),size=guide_legend(ncol=2))+
  #scale_color_gradient(low = "#0766AD",high = "#6EC207",guide = "colorbar")+
  scale_color_gradient(low = "#1f96c4",high = "#20684a",guide = "colorbar")+
  labs(x=paste0("PCoA 1 (",round((pcoa_occ$eig[1]/sum(pcoa_occ$eig))*100,1),"%)"),y=paste0("PCoA 2 (",round((pcoa_occ$eig[2]/sum(pcoa_occ$eig))*100,1),"%)"),color="Proportion\nof terrestrial\ngroups",size="Proportion\nof forest")+
  theme(legend.justification = "left",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=9), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));CompPlt2

##### Taxon-specific (Abundance) ----
taxa_mat = melt(ins_mat)
colnames(taxa_mat) = c("Trap","Taxa","Abun")

gam_data<-data.frame(taxa_mat,Long=land$X,Lat=land$Y,Area=land[,"Area0100"])

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
  
  return(list(moran=data.frame(observed=MoranI$observed,expected=MoranI$expected,sd=MoranI$sd,p=MoranI$p.value, taxa=w),
              coefs=data.frame(Estimate=c(summary(GAM_Abun)[[1]],NA),k=temp_k,stat=c(summary(GAM_Abun)[[3]],summary(GAM_Abun)[[7]]),P=c(summary(GAM_Abun)[[4]],summary(GAM_Abun)[[8]]),R2=summary(GAM_Abun)[[10]],AIC=AIC(GAM_Abun),var=c(names(summary(GAM_Abun)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_Abun)[[1]])),"Random effects"), taxa=w)))
})
do.call("rbind",lapply(taxa_gam_res,"[[",1))->moran_taxa
do.call("rbind",lapply(taxa_gam_res,"[[",2))->SoE_taxa

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

taxa_mat$Area<-land$Area0100

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

plt_data = taxa_mat[taxa_mat$Abun<500,]
plt_data1 = plt_data[grepl("Caddisflies|Mayflies|Mosquitoes",plt_data$Taxa),]
plt_data2 = plt_data[!grepl("Caddisflies|Mayflies|Mosquitoes",plt_data$Taxa),]

AbunOrdPlt1<-ggplot()+
  geom_point(data=plt_data,aes(x=Area,y=log10(Abun+1),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=plt_data,aes(x=Area,y=log10(Abun+1),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
  facet_wrap(.~Taxa,scales = "free_y",nrow = 2)+
  scale_color_manual(values=my_cols,labels=fantasy_names)+
  scale_fill_manual(values=my_cols,labels=fantasy_names)+
  scale_linetype(labels=fantasy_names)+
  theme_minimal()+
  scale_y_continuous(labels = c(0,10,100),breaks = c(0,1,2))+
  #guides(color="none",fill="none",color="none",shape="none")+
  scale_shape_manual(values=c(21,22,23,24,25,8,13),labels=fantasy_names)+
  labs(x="Proportion of forest",y="Number of individuals",color="Biological\ngroup",fill="Biological\ngroup",shape="Biological\ngroup",linetype="Biological\ngroup")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));AbunOrdPlt1

#### Taxon-specific (Body size) ----
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
gam_data<-data.frame(size_mat,Long=land$X,Lat=land$Y,Area=land[,"Area0100"])

taxa_gam = gam_data
taxa_gam = na.omit(taxa_gam)
levels(taxa_gam$Taxa) = as.character(unique(taxa_gam$Taxa))[order(as.character(unique(taxa_gam$Taxa)))]
AICs<-lapply(1:20,function(y){
  gam(Size~Area*Taxa+s(Long, Lat,k=y),family = "quasipoisson",data = taxa_gam)->GAM_Size
  return(AIC(GAM_Size))
})
data.frame(k=1:20,AIC=unlist(AICs))->k_data
k_data[order(k_data$AIC),][1,]$k->temp_k
glm(Size~Area+Taxa+Area:Taxa,family = "quasipoisson",data = taxa_gam)->GAM_Size
summary(GAM_Size)
residuals(GAM_Size)->size_res
size_ins = aggregate(size_res,list(taxa_gam$Taxa),mean)[,2]
names(size_ins) = unique(insects$class)[order(unique(insects$class))]
K_result <- phylosig(tree.ins, size_ins, method = "lambda", test = TRUE,niter = 9999)
print(K_result)

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
  
  return(list(moran=data.frame(observed=MoranI$observed,expected=MoranI$expected,sd=MoranI$sd,p=MoranI$p.value, taxa=w),
              coefs=data.frame(Estimate=c(summary(GAM_Size)[[1]],NA),k=temp_k,stat=c(summary(GAM_Size)[[3]],summary(GAM_Size)[[7]]),P=c(summary(GAM_Size)[[4]],summary(GAM_Size)[[8]]),R2=summary(GAM_Size)[[10]],AIC=AIC(GAM_Size),var=c(names(summary(GAM_Size)[[1]]),"Coordinates"),eff=c(rep("Fixed effects",length(summary(GAM_Size)[[1]])),"Random effects"), taxa=w)))
})
do.call("rbind",lapply(taxa_gam_res,"[[",1))->moran_taxa
do.call("rbind",lapply(taxa_gam_res,"[[",2))->SoE_size
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

size_mat$Area<-land$Area0100
  
as.character(size_mat$Taxa)->size_mat$Taxa
size_mat$Taxa = factor(size_mat$Taxa,levels=all_names$Taxa)

as.character(area_coefs_size_ord$taxa)->area_coefs_size_ord$taxa
area_coefs_size_ord$taxa = factor(area_coefs_size_ord$taxa,levels=all_names$Taxa)

as.character(area_size$taxa)->area_size$taxa
area_size$taxa = factor(area_size$taxa,levels=all_names$Taxa)

size_mat = size_mat[!is.na(size_mat$Size),]
size_mat[(size_mat$Taxa!="Trichoptera")&(size_mat$Taxa!="Hemiptera")&(size_mat$Taxa!="Nematocera"),]->size_plt
size_plt$Taxa<-all_names[match(size_plt$Taxa,all_names$Taxa),2]
size_plt$Taxa = factor(size_plt$Taxa,levels = c("Bees and wasps","Beetles","Flies","Mayflies"))

size_plt$Taxa = paste0("(",letters[1:5],"): ",levels(size_plt$Taxa))[match(size_plt$Taxa,levels(size_plt$Taxa))]

SizeOrdPlt1<-ggplot()+
  geom_point(data=size_plt,aes(x=Area,y=log10(Size),color=Taxa,fill=Taxa,shape=Taxa))+
  geom_smooth(data=size_plt,aes(x=Area,y=log10(Size),color=Taxa),method = "glm",method.args=list(family="quasipoisson"))+
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

##### Arrange and save plots ----
ggsave(filename = "figures/abundance_order.tif",plot = AbunOrdPlt1,units = "in",width = 13.5,height = 6,dpi = 600)

ggsave(filename = "figures/size_order.tif",plot = SizeOrdPlt1,units = "in",width = 7,height = 6,dpi = 600)

arr3 = ggarrange(SADPlt,NPlt,common.legend = T,legend = "right");arr3
ggsave(filename = "figures/SAD+N.tif",plot = arr3,units = "in",width = 9.5,height = 4,dpi = 600)

arr6 = ggarrange(betaC_plt,ggarrange(CompPlt2,CompPlt1,common.legend = T,legend = "right"),nrow = 1,ncol = 2,widths = c(0.275,1),align = "h");arr6
ggsave(filename = "figures/composition_plot.tif",plot = arr6,units = "in",width = 15,height = 5,dpi = 600)

##### Supporting figures ----
dir("model_training/runs/",pattern = "^train",full.names = T) -> mod_files
lapply(paste0(mod_files,"/results.csv"),read.csv,h=T)->mod_res
for(x in 1:length(mod_res)){
  mod_res[[x]]$fold<-x
}
do.call("rbind",mod_res)->mod_res
melt(mod_res,id.vars = c("epoch","fold"))->mod_res
mod_res = mod_res[mod_res$variable!="time",]
mod_res = mod_res[mod_res$variable!="lr.pg1",]
mod_res = mod_res[mod_res$variable!="lr.pg2",]
facet_labs=c("Train: box loss","Train: class. loss","Train: dfl. loss","Train: box precision","Train: recall","Train: mAP (@0.5 IoU)","Train: mAP (@0.95 IoU)","Val.: box loss","Val.: class. loss","Val.: dfl. loss","Learning rate")
names(facet_labs) = unique(mod_res$variable) 
metrics_plt = ggplot(data = mod_res,mapping = aes(x=epoch,y=value))+
  geom_line(aes(group=factor(fold),color=factor(fold),alpha=factor(fold)))+
  geom_smooth(color="black",se=F,method = "loess")+
  facet_wrap(.~variable,scales = "free_y",labeller = labeller(variable = facet_labs))+
  theme_minimal()+
  scale_x_continuous(limits = c(0,250))+
  scale_alpha_manual(values=c(rep(0.5,5)),labels=c("1","2","3","4","5"),name="Fold")+
  scale_color_manual(values=c("#BB3E00","#F7AD45","#657C6A","#A2B9A7","#547792"),labels=c("1","2","3","4","5"),name="Fold")+
  labs(y="Value",x="Epoch")+
  theme(legend.justification = "top",strip.text = element_text(size=12,family="sans",face = "bold",hjust = 0,vjust = 0.5),strip.background = element_rect(fill="transparent",color="transparent"), legend.text.align = 1,legend.title.align=1,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));metrics_plt
ggsave(filename = "figures/metrics_plot.tif",plot = metrics_plt,units = "in",width = 12.5,height = 7.5,dpi = 600)

aggregate(mod_res[mod_res$epoch==250,]$value,list(mod_res[mod_res$epoch==250,]$variable),mean)->mod_res_agg
aggregate(mod_res[mod_res$epoch==250,]$value,list(mod_res[mod_res$epoch==250,]$variable),sd)->mod_res_agg_sd

### Arrange Moran's Table

TableS2 = data.frame(Group=c(rep("Lifecycle",nrow(moran_S)+nrow(moran_Comp)),as.character(moran_taxa$taxa),as.character(moran_size$taxa)),
           Observed=c(moran_S$observed,moran_Comp$observed,moran_taxa$observed,moran_size$observed),
           Expected=c(moran_S$expected,moran_Comp$expected,moran_taxa$expected,moran_size$expected),
           SD=c(moran_S$sd,moran_Comp$sd,moran_taxa$sd,moran_size$sd),
           P=c(moran_S$p,moran_Comp$p,moran_taxa$p,moran_size$p),
           Scale=c(moran_S$scale,moran_Comp$scale,moran_taxa$scale,moran_size$scale),
           Response=c(moran_S$component,moran_Comp$data,rep("Abundance",nrow(moran_taxa)),rep("Body size",nrow(moran_size))))
TableS2[order(TableS2$Response),]->TableS2
write.csv(TableS2,"results/TableS2.csv")

gsub("habT","Life cycle",SoE_S$var)->SoE_S$var
TableS3 = data.frame(Response=c(paste0(SoE_S$component,"-component"),paste0(SoE_Comp$data,"-based composition")),
                     "Effect type"=c(SoE_S$eff,SoE_Comp$eff),
                     Explanatory=c(SoE_S$var,SoE_Comp$var),
                     Estimate=c(SoE_S$Estimate,SoE_Comp$Estimate),
                     Statistics=c(SoE_S$stat,SoE_Comp$stat),
                     K=c(SoE_S$k,SoE_Comp$k),
                     P=c(SoE_S$P,SoE_Comp$P),
                     Signif=c(SoE_S$P.signif,SoE_Comp$P.signif),
                     R2=c(SoE_S$R2,SoE_Comp$R2),
                     AIC=c(SoE_S$AIC,SoE_Comp$AIC),
                     Scale=c(SoE_S$scale,SoE_Comp$scale))
TableS3
write.csv(TableS3,"results/TableS3.csv")

TableS4 = data.frame(Response=c(rep("Abundance",nrow(SoE_taxa)),rep("Body size",nrow(SoE_size))),
                     Taxon=c(SoE_taxa$taxa,SoE_size$taxa),
                     "Effect type"=c(SoE_taxa$eff,SoE_size$eff),
                     Explanatory=c(SoE_taxa$var,SoE_size$var),
                     Estimate=c(SoE_taxa$Estimate,SoE_size$Estimate),
                     Statistics=c(SoE_taxa$stat,SoE_size$stat),
                     K=c(SoE_taxa$k,SoE_size$k),
                     P=c(SoE_taxa$P,SoE_size$P),
                     Signif=c(SoE_taxa$P.signif,SoE_size$P.signif),
                     R2=c(SoE_taxa$R2,SoE_size$R2),
                     AIC=c(SoE_taxa$AIC,SoE_size$AIC),
                     Scale=c(SoE_taxa$scale,SoE_size$scale))
TableS4
write.csv(TableS4,"results/TableS4.csv")
