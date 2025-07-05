# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # Data extraction script for Colares et al., 2025 # # # # #
# # # # # # # # # # # # # # By L. Colares # # # # # # # # # # # # #
source("scripts/00. setup.R")

##### 1. Biological data ----
##### 1.1. Load detectiosn from YOLOv8 ----
read.table("datasets/DL-fold/1-fold/data.yaml",quote = "",fill = T)->label_map
data.frame(ID=seq(0,length(label_map[5,2:ncol(label_map)])-1),group=gsub("\\[|'|\\]|,","",unlist(label_map[5,2:ncol(label_map)])))->label_map
rownames(label_map)=1:nrow(label_map)

detected = dir(paste0("datasets/IDs/",1:5,"-fold"),pattern = "*.csv$",full.names = T)
detc = lapply(1:length(detected),function(x){
  read.csv(detected[x])[,-1]->detc
  detc$img<-tail(strsplit(detected[x],"/")[[1]],1)
  detc$fold<-strsplit(detected[x],"/")[[1]][grepl("fold",strsplit(detected[x],"/")[[1]])]
  return(detc)
})
do.call("rbind",detc)->detc
detc$fold = as.numeric(gsub("-fold","",detc$fold))
detc$img = gsub("\\.csv","",detc$img)
detc$trap = substr(detc$img,1,nchar(detc$img)-2)

##### 1.1. Calculate confidence thresholds for each group ----
val_det = dir(paste0("model_training/runs/val",1:5,"/labels"),pattern = "*.txt$",full.names = T)
load_val = lapply(1:length(val_det),function(x){
  read.table(val_det[x])->det
  det$img<-tail(strsplit(val_det[x],"/")[[1]],1)
  det$fold<-strsplit(val_det[x],"/")[[1]][grepl("val",strsplit(val_det[x],"/")[[1]])]
  return(det)
})
do.call("rbind",load_val)->load_val
load_val$fold = as.numeric(gsub("val","",load_val$fold))
load_val$img = gsub("\\.txt","",load_val$img)

true_val = dir(paste0("datasets/DL-fold/",1:5,"-fold/val/labels"),pattern = "*.txt$",full.names = T)
true_val = lapply(1:length(true_val),function(x){
  read.table(true_val[x])->det2
  det2$img<-tail(strsplit(true_val[x],"/")[[1]],1)
  det2$fold<-strsplit(true_val[x],"/")[[1]][grepl("^[0-9]-fold",strsplit(true_val[x],"/")[[1]])]
  return(det2)
})
do.call("rbind",true_val)->true_val
true_val$fold = as.numeric(gsub("-fold","",true_val$fold))
true_val$img = gsub("\\.txt","",true_val$img)

data.frame(matrix(NA,nrow = length(seq(0,1,0.025)),ncol = 7))->metrics
colnames(metrics)=c("Thresh","TP","FP","FN","Accurracy","Recall","F1")
metrics$Thresh<-seq(0,1,0.025)

load_val_copy=load_val

val_metric_fold<-lapply(1:5,function(fold){
  val_metric<-pblapply(unique(true_val$img),function(img.name){
    if((nrow(load_val[(load_val$img==img.name)&(load_val$fold==fold),])!=0)&(nrow(true_val[(true_val$img==img.name)&(true_val$fold==fold),])!=0)){
      dupl_yolo=rbind(data.frame(true_val[(true_val$img==img.name)&(true_val$fold==fold),1:5],V6=NA,cat="true"), data.frame(load_val[(load_val$img==img.name)&(load_val$fold==fold),1:6],cat="val"))
      dupl_yolo[,2:5] = dupl_yolo[,2:5]*640
      
      yolo_format = matrix(NA,nrow = nrow(dupl_yolo),ncol = 4)
      colnames(yolo_format)=c("xmin","ymin","xmax","ymax")
      yolo_format[,1]=dupl_yolo$V2-(dupl_yolo$V4/2)
      yolo_format[,2]=dupl_yolo$V3-(dupl_yolo$V5/2)
      yolo_format[,3]=dupl_yolo$V2+(dupl_yolo$V4/2)
      yolo_format[,4]=dupl_yolo$V3+(dupl_yolo$V5/2)
      dupl_yolo[,2:5]=yolo_format
      
      colnames(dupl_yolo)=c("class","xmin","ymin","xmax","ymax","conf","cat")
      
      dupl_val = dupl_yolo[dupl_yolo$cat=="val",]
      rownames(dupl_val)=1:nrow(dupl_val)
      
      IoU_val<-IoU(dupl_val)
      IoU_val=as.matrix(IoU_val)
      rownames(IoU_val)=1:nrow(IoU_val)
      pair_comp<-reshape2::melt(IoU_val)
      
      non_over_temp<-{}
      
      for(u in 1:nrow(dupl_val)){
        sel_fis<-pair_comp[pair_comp$Var1==rownames(dupl_val[u,]),]
        sel_fis[sel_fis$value>=0.5,]->sel_fis
        if(nrow(sel_fis)>0){
          over_bbox<-unique(c(sel_fis$Var1,sel_fis$Var2))
          non_over_temp[[length(non_over_temp)+1]]<-dupl_val[grepl(paste0("^",over_bbox,"$",collapse = "|"),rownames(dupl_val)),][order(dupl_val[grepl(paste0("^",over_bbox,"$",collapse = "|"),rownames(dupl_val)),]$conf,decreasing = T),][1,]
        } else {
          dupl_val[u,]->non_over_temp[[length(non_over_temp)+1]]
        } 
      }
      do.call("rbind",non_over_temp)->non_over_temp
      non_over_temp[!duplicated(non_over_temp[,1:5]),]->non_over_bbox
      
      dupl_yolo[dupl_yolo$cat=="true",]->true_data
      IoU4Acc<-suppressMessages(IoU2(mat1 = true_data, mat2 = non_over_bbox))
      colnames(IoU4Acc)=c("true","val","IoU")
      IoU4Acc<-as.data.frame(IoU4Acc)
      
      non_over_bbox$result<-NA
      true_data$result<-NA
      
      true_data[match(IoU4Acc[IoU4Acc$IoU>0.5,]$true,rownames(true_data)),]$result<-true_data[match(IoU4Acc[IoU4Acc$IoU>0.5,]$true,rownames(true_data)),]$class==non_over_bbox[match(IoU4Acc[IoU4Acc$IoU>0.5,]$val,rownames(non_over_bbox)),]$class
      
      true_data$pred_lab<-NA
      true_data[match(IoU4Acc[IoU4Acc$IoU>0.5,]$true,rownames(true_data)),]$pred_lab<-non_over_bbox[match(IoU4Acc[IoU4Acc$IoU>0.5,]$val,rownames(non_over_bbox)),]$class
      true_data[match(IoU4Acc[IoU4Acc$IoU>0.5,]$true,rownames(true_data)),]$conf<-non_over_bbox[match(IoU4Acc[IoU4Acc$IoU>0.5,]$val,rownames(non_over_bbox)),]$conf
      
      new_non_over<-non_over_bbox[-match(IoU4Acc[IoU4Acc$IoU>0.5,]$val,rownames(non_over_bbox)),]
      
      if(nrow(new_non_over)!=0){
        new_non_over$pred_lab<-NA
        final_met<-rbind(true_data,new_non_over)
        
        
        if(nrow(final_met[is.na(final_met$result)&final_met$cat=="val",])!=0){
          final_met[is.na(final_met$result)&final_met$cat=="val",]$result<-"FP"
        }
        if(nrow(final_met[is.na(final_met$result)&final_met$cat=="true",])!=0){
          final_met[is.na(final_met$result)&final_met$cat=="true",]$result<-"FN"
        }
        if(nrow(final_met[final_met$result=="TRUE",])!=0){
          final_met[final_met$result=="TRUE",]$result<-"TP"
        }
        if(nrow(final_met[final_met$result=="FALSE",])!=0){
          final_met[final_met$result=="FALSE",]$result<-"FP"
        }
        
        return(data.frame(final_met,img=img.name,fold=fold))   
      } else {
        final_met<-true_data
        
        
        if(nrow(final_met[is.na(final_met$result)&final_met$cat=="val",])!=0){
          final_met[is.na(final_met$result)&final_met$cat=="val",]$result<-"FP"
        }
        if(nrow(final_met[is.na(final_met$result)&final_met$cat=="true",])!=0){
          final_met[is.na(final_met$result)&final_met$cat=="true",]$result<-"FN"
        }
        if(nrow(final_met[final_met$result=="TRUE",])!=0){
          final_met[final_met$result=="TRUE",]$result<-"TP"
        }
        if(nrow(final_met[final_met$result=="FALSE",])!=0){
          final_met[final_met$result=="FALSE",]$result<-"FP"
        }
        
        return(data.frame(final_met,img=img.name,fold=fold))   
      }
    }
  })
  return(do.call("rbind",Filter(Negate(is.null), val_metric)))
})
do.call("rbind",val_metric_fold)->val_metric_fold

val_metric_fold[is.na(val_metric_fold$conf),]$conf<-1

all_metrics_fold<-lapply(1:5,function(fold){
  all_metrics<-pblapply(c(seq(0,0.99,0.01)),function(thresh){
    #print(thresh)
    val_metric_copy=val_metric_fold[val_metric_fold$fold==fold,]
    if(sum((val_metric_copy$conf<thresh)&(val_metric_copy$result=="TP"))>0){
      val_metric_copy[(val_metric_copy$conf<thresh)&(val_metric_copy$result=="TP"),]$result<-"FN"
    }
    if(sum((val_metric_copy$conf<thresh)&(val_metric_copy$result=="FP"))>0){
      val_metric_copy[(val_metric_copy$conf<thresh)&(val_metric_copy$result=="FP"),]$result<-"TN"
    }
    
    mets<-table(val_metric_copy$result)
    
    if(sum(grepl("TP",names(mets)))==0){
      mets<-c(mets,0)
      names(mets)=gsub("^$","TP",names(mets))
    }
    if(sum(grepl("FN",names(mets)))==0){
      mets<-c(mets,0)
      names(mets)=gsub("^$","FN",names(mets))
    }
    if(sum(grepl("FP",names(mets)))==0){
      mets<-c(mets,0)
      names(mets)=gsub("^$","FP",names(mets))
    }
    
    precision<-mets["TP"]/(mets["TP"]+mets["FP"])
    recall<-mets["TP"]/(mets["TP"]+mets["FN"])
    F1<-(2*precision*recall)/(precision+recall)
    f1_class<-lapply(unique(val_metric_copy$class),function(y){
      mets2<-table(val_metric_copy[val_metric_copy$class==y,]$result)
      if(sum(grepl("TP",names(mets2)))==0){
        mets2<-c(mets2,0)
        names(mets2)=gsub("^$","TP",names(mets2))
      }
      if(sum(grepl("FN",names(mets2)))==0){
        mets2<-c(mets2,0)
        names(mets2)=gsub("^$","FN",names(mets2))
      }
      if(sum(grepl("FP",names(mets2)))==0){
        mets2<-c(mets2,0)
        names(mets2)=gsub("^$","FP",names(mets2))
      }
      precision2<-mets2["TP"]/(mets2["TP"]+mets2["FP"])
      recall2<-mets2["TP"]/(mets2["TP"]+mets2["FN"])
      F12<-(2*precision2*recall2)/(precision2+recall2)
      return(data.frame(Thresh=thresh,precision=precision2,recall=recall2,F1=F12,group=y))
    })
    return(rbind(data.frame(Thresh=thresh,precision=precision,recall=recall,F1=F1,group="all"),do.call("rbind",f1_class)))
  })
  do.call("rbind",all_metrics)->all_metrics
  rownames(all_metrics)=1:nrow(all_metrics)
  all_metrics[is.na(all_metrics)]<-0
  return(data.frame(all_metrics,fold=fold))
})
do.call("rbind",all_metrics_fold)->all_metrics_fold

library(ggplot2)
all_metrics_fold$group<-label_map[match(all_metrics_fold$group,label_map$ID),]$group
all_metrics_fold$group[is.na(all_metrics_fold$group)]<-"All groups"
ggplot(data = all_metrics_fold[all_metrics_fold$group!="Orthoptera",],mapping = aes(x=Thresh,y=F1,color=group,group=group))+
  geom_line()+
  facet_wrap(.~fold)

bestThresh<-lapply(unique(all_metrics_fold$fold), function(x){
  group_best<-lapply(unique(all_metrics_fold$group), function(y){
    sels<-all_metrics_fold[(all_metrics_fold$fold==x)&(all_metrics_fold$group==y),]
    return(sels[order(sels$precision,decreasing = T),][1,])
  })
  return(do.call("rbind",group_best))
})
do.call("rbind",bestThresh)->bestThresh
bestThresh = bestThresh[bestThresh$group!="All groups",]

bestThreshF1<-lapply(unique(all_metrics_fold$fold), function(x){
  group_best<-lapply(unique(all_metrics_fold$group), function(y){
    sels<-all_metrics_fold[(all_metrics_fold$fold==x)&(all_metrics_fold$group==y),]
    return(sels[order(sels$F1,decreasing = T),][1,])
  })
  return(do.call("rbind",group_best))
})
do.call("rbind",bestThreshF1)->bestThreshF1
write.csv(bestThreshF1,"results/metrics_at_F1.csv")
mean(bestThreshF1[bestThreshF1$group=="All groups",]$precision)

##### 1.2. Apply best threshold to each group and fold ----
final_mat = {}
for (x in unique(detc$fold)) {
  for(y in unique(detc$class)){
    final_mat[[length(final_mat)+1]] = detc[(detc$fold==x)&(detc$class==y)&(detc$conf>=bestThresh[(bestThresh$fold==x)&(bestThresh$group==y),]$Thresh),]
  }
}
final_mat = do.call("rbind",final_mat)

true_val$V1 = label_map[match(true_val$V1,label_map$ID),2]

unique(c("Araneae",bestThresh[bestThresh$Thresh==0,]$group,names(table(true_val$V1))[table(true_val$V1)<=30]))->rem_group

final_mat = final_mat[!grepl(paste0("^",rem_group,"$",collapse = "|"),final_mat$class),]
rownames(final_mat)=1:nrow(final_mat)

final_ind = pblapply(unique(final_mat$img), function(y){
  print(y)
  final_mat[final_mat$img==y,]->img_mat
  rownames(img_mat) = 1:nrow(img_mat)
  
  overlaps<-as.matrix(IoU(img_mat))
  
  reshape2::melt(overlaps)->overlaps
  overlaps[overlaps$Var1!=overlaps$Var2,]->overlaps
  overlaps[overlaps$value>=0.5,]->overlaps
  
  if(nrow(overlaps)!=0){
    rem_rows<-lapply(unique(overlaps$Var1),function(x){
      unique(unlist(overlaps[(overlaps$Var1==x)|(overlaps$Var2==x),1:2]))->sel_row
      return(as.numeric(rownames(img_mat[sel_row,][order(img_mat[sel_row,]$conf,decreasing = T),][2:nrow(img_mat[sel_row,][order(img_mat[sel_row,]$conf,decreasing = T),]),])))
    })
    
    unlist(rem_rows)->rem_rows
    unique(rem_rows)->rem_rows
    img_mat = img_mat[-rem_rows,]
    return(img_mat)  
  } else {
    return(img_mat)
  }
  
})
final_ind2 = do.call("rbind",final_ind)
write.csv(final_ind2,"matrix of insects.csv")

##### 2. Extract landscape metrics -----
data.frame(read.csv("datasets/csv/coordinates.csv",h=T,sep = ";"))->coords
coords = st_as_sf(coords,coords = c("Long","Lat"),crs = 4326)
coords_t = st_transform(coords,crs = "+proj=utm +zone=21 +south +ellps=GRS80 +units=m +no_defs")

download.file("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_9/lclu/coverage/brasil_coverage_2021.tif",destfile = "datasets/spatial/brasil_coverage_2021.tif")
c(seq(100,1000,50),seq(1100,10000,100))->buf_sizes

raster("datasets/spatial/brasil_coverage_2021.tif")->mapbiomas
crop(mapbiomas,st_buffer(coords,5000))->map_cropped
rm(mapbiomas)
projectRaster(map_cropped,crs = "+proj=utm +zone=21 +south +ellps=GRS80 +units=m +no_defs",method = "ngb")->map_cropped
lsm_c_ca(map_cropped)->all_classes
all_classes<-data.frame(class=all_classes$class,cat=c("NA","land","land","land","land","land","water","land"))

snow::makeCluster(4)->cl
snow::clusterExport(cl,list(c("map_cropped","all_classes")))
final_area<-lapply(buf_sizes,function(x){
  st_buffer(coords_t,x)->bufs
  snow::clusterExport(cl,list("bufs"))
  print(paste0(x," meters:"))
  areaz<-pblapply(cl = cl,1:nrow(bufs),function(y){
    raster::crop(map_cropped,raster::extent(bufs[y,]))->map_buf
    raster::mask(map_buf,bufs[y,])->map_buf
    landscapemetrics::lsm_c_ca(map_buf)->class_area
    class_area$value<-class_area$value/sum(class_area$value)
    return(sum(class_area[na.omit(match(all_classes$class[all_classes$cat=="land"],class_area$class)),]$value))
  })
  area_data<-data.frame(Area=unlist(areaz))
  colnames(area_data)=sprintf("%04d",x)
  return(area_data)
})
do.call("cbind",final_area)->final_area
colnames(final_area)=paste0("Area",colnames(final_area))
final_area = data.frame(Trap=coords$Trap,st_coordinates(coords),final_area)
write.csv(final_area,"datasets/csv/landscape metrics - jul25.csv")
