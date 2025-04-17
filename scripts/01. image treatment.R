# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # Image treatment script for Colares et al., 2025 # # # # #
# # # # # # # # # # # # # # By L. Colares # # # # # # # # # # # # #

##### 00. Download raw images from FigShare ----

download.file(url = "https://figshare.com/ndownloader/files/41804211",
              destfile = "datasets/imgs/originals.zip",mode = "wb")
unzip(zipfile = "datasets/imgs/originals.zip",exdir = "datasets/imgs/")
file.rename("datasets/imgs/originals","datasets/imgs/raw-imgs")

# If successful, you should have three folders with raw, sliced and labelled, which are the same as the ones used in the published paper. You can skip to step 03 if this is the case.

##### 01. Slice images for annotation ----

images=slices=dir("datasets/imgs/raw-imgs/",full.names = T)
destination.folder="datasets/imgs/new-slices/"
dir.create(destination.folder)
mask_segmentation(imgs = images[1],destFolder = destination.folder,Resize = 4, K=0.2, Blur=5, H=0.05, Margin=100)

##### 02. Pre-annotate with the interactive function ----
## Alternatively, use labelImg to annotate the images or go to the next step if you already have the data

slices=dir("datasets/imgs/sliced-imgs/",full.names = T) # or slices=dir("datasets/imgs/new-slices/",full.names = T)
destination.folder2="datasets/imgs/new-annotated/"
dir.create(destination.folder2)
auto_annotation(imgs = slices, destFolder = destination.folder2, Blur = 2, K = 0.2, Erode = 5)

##### 03. Split dataset into train and validation across 5 folds ----

pascal_ann = dir("datasets/imgs/labelled-imgs/",pattern = "*.xml$",full.names = T)
pascal_img = dir("datasets/imgs/labelled-imgs/",pattern = "*.jpg$",full.names = T)

dir.create("datasets/DL-fold/base-fold")
yolo_anns = pascal_to_yolo(ann = pascal_ann,img = pascal_img,resize = T,size = 640,flat = T,save_ann = T,save_img = T,savePath = paste0("datasets/DL-fold/base-fold"))

yolo_ann = dir("datasets/DL-fold/base-fold/",pattern = "*.txt$",full.names = T)
yolo_img = dir("datasets/DL-fold/base-fold/",pattern = "*.jpg$",full.names = T)
data.frame(ann=yolo_ann,img=yolo_img,priority=NA,fold=NA)->final_fold

all_anns = lapply(1:length(yolo_ann),function(x){
  load_ann = read.table(yolo_ann[x])
  return(data.frame(load_ann,img=yolo_img[x]))
})
do.call("rbind",all_anns)->all_anns
table(all_anns$V1)
data.frame(names=names(table(all_anns$V1)[order(table(all_anns$V1))]),priority=1:length(table(all_anns$V1)))->priority_set
all_anns$fold<-NA
all_anns$priority<-NA

for(x in unique(all_anns$img)){
  priority_set[match(all_anns[all_anns$img==x,]$V1,priority_set$names),][order(priority_set[match(all_anns[all_anns$img==x,]$V1,priority_set$names),]$priority),1][1]->priority_group
  all_anns[all_anns$img==x,]$priority<-as.numeric(priority_group)
}

all_anns[match(final_fold$img,all_anns$img),]$priority->final_fold$priority

splitted<-lapply(unique(final_fold$priority), function(x){
  final_fold[as.numeric(sample(rownames(final_fold[final_fold$priority==x,]))),]->sampled_fold
  split(sampled_fold, rep(1:5,each=1))->splitted
  for(y in 1:length(splitted)){
    splitted[[y]]$fold<-y
  }
  return(do.call("rbind",splitted))
})
do.call("rbind",splitted)->splitted

savePath = "datasets/DL-fold/"
for(x in 1:5){
  dir.create(file.path(savePath,paste0(x,"-fold")))
  
  dir.create(file.path(savePath,paste0(x,"-fold"),"train"))
  dir.create(file.path(savePath,paste0(x,"-fold"),"train","images"))
  dir.create(file.path(savePath,paste0(x,"-fold"),"train","labels"))
  
  dir.create(file.path(savePath,paste0(x,"-fold"),"val"))
  dir.create(file.path(savePath,paste0(x,"-fold"),"val","images"))
  dir.create(file.path(savePath,paste0(x,"-fold"),"val","labels"))
  
  file.copy(splitted[splitted$fold!=x,]$ann,file.path(savePath,paste0(x,"-fold"),"train","labels",unlist(lapply(strsplit(splitted[splitted$fold!=x,]$ann,"/"),tail,1))))
  file.copy(splitted[splitted$fold!=x,]$img,file.path(savePath,paste0(x,"-fold"),"train","images",unlist(lapply(strsplit(splitted[splitted$fold!=x,]$img,"/"),tail,1))))
  
  file.copy(splitted[splitted$fold==x,]$ann,file.path(savePath,paste0(x,"-fold"),"val","labels",unlist(lapply(strsplit(splitted[splitted$fold==x,]$ann,"/"),tail,1))))
  file.copy(splitted[splitted$fold==x,]$img,file.path(savePath,paste0(x,"-fold"),"val","images",unlist(lapply(strsplit(splitted[splitted$fold==x,]$img,"/"),tail,1))))
  
  YOUtils::save_label_map(labmap = yolo_anns$LabelMap,savePath = file.path(savePath,paste0(x,"-fold")))
}

#### Next steps: train the YOLOv8 model using the provided Python scripts, either locally or using Google Colab
# End of script