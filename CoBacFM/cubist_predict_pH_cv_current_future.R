OTUcoreS1 <- read.table("/data/fengkai/grassland/UNOISE/CoreSpecies/S1_core_otu_relative_abundance_989_species_1251_samples.txt",sep="\t",row.names = 1,header=T)
RawEnv <-read.table("/data/fengkai/grassland/UNOISE/CoreSpecies/1251samples_env_table.txt",sep="\t",header=TRUE,row.names=1)

SampleUniq <- intersect(colnames(OTUcoreS1),rownames(RawEnv))
OTUcoreS1 <- OTUcoreS1[,SampleUniq]
ENV1251 <- RawEnv[SampleUniq,]

sle <- c("pH","SOC","TN","BIO1","BIO5","BIO12","BIO13","plant_ndvi_2010_2020")
library(vegan)
env.st <- ENV1251[,sle]
st <- FALSE # TRUE,FALSE
if(st){
  env.sd <- decostand(env.st,method = "standardize",na.rm=T)
}else{
  env.sd <- env.st
}
colnames(env.sd) <- c("pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")

eco.tab <- read.csv("/data/fengkai/grassland/UNOISE/RF/selected_reduced_699 OTUs cluster names.csv",sep=",",header = T,row.names = 1)
eco.grp <- list()
for (i in 1:length(unique(eco.tab[,ncol(eco.tab)]))) {
  eco.grp[[i]] <- rownames(eco.tab)[which(eco.tab[,ncol(eco.tab)]==unique(eco.tab[,ncol(eco.tab)])[i])]
}
names(eco.grp) <- unique(eco.tab[,ncol(eco.tab)])
eco <- names(eco.grp) 
eco <- eco[-c(which(eco==""))]

eco.abund <- c()
for(slc in eco){
  eco.abund <- cbind(eco.abund,colSums(OTUcoreS1[eco.grp[[slc]],]))
}
colnames(eco.abund) <- eco
sum(rownames(eco.abund)==rownames(env.sd))

reg.method <-  c("cubist")
standard <- c("Relative")
for(fm in c("BCC_CSM2_MR","IPSL_CM6A_LR")){
  for(sn in c('ssp126','ssp370','ssp585')){
    for(y in c("2021","2041","2061","2081")){
      eco.future <- c()
      for (slc in eco) {
        cluster.fu <- read.table(paste("/data/fengkai/grassland/UNOISE/future/future_predicted",fm,sn,y,standard, "abundance",slc,reg.method,"value.txt",sep="_"),
                                 header=T,row.names = 1,sep="\t")
        cluster.fu.sl <- cluster.fu[,c("long","lat","pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI","mean")]
        geo <- cluster.fu.sl[,1:10]
        eco.future <- cbind(eco.future,cluster.fu.sl[,11])
      }
      eco.future <- cbind(geo,eco.future)
      colnames(eco.future) <- c("long","lat","pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI",eco)
      write.table(eco.future,paste("/data/fengkai/grassland/UNOISE/future/Clusters_future_predicted",fm,sn,y,standard, "abundance",reg.method,"value.txt",sep="_"),sep="\t",col.names=NA,row.names = T,quote=F)
    }
  }
}

env_GC30_grass_current <- read.table("/data/fengkai/gis/GlobalLand30/GC30Env_current_extracted.txt",header=T,row.names = 1,sep="\t")
if(!dir.exists("/data/fengkai/grassland/UNOISE/future/pH")){dir.create("/data/fengkai/grassland/UNOISE/future/pH")}
standard = "Relative"
cvnum = 10 #　10-fold cross validation
for(reg.method in c("cubist")){ # "randomforest",
  pre.dat <- cbind(eco.abund,env.sd[,c("pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")]) 
  if(standard == "Standardized"){
    pre.dat[,ncol(pre.dat)] <- decostand(pre.dat[,ncol(pre.dat)],method = "standardize") 
  }
  pre.dat <- pre.dat[,c("Low pH","High pH and TN with low SOC","Low SOC","High SOC","High PWetM","Low pH high SOC","High MaxTWarmM",
                        "Low pH low MAT","High pH","High NDVI","pH")]  # select pH as predicting variable
  set.seed(36)
  neworder <- sample(1:nrow(pre.dat),nrow(pre.dat),replace = FALSE)
   for(fm in c("BCC_CSM2_MR","IPSL_CM6A_LR")){
    for(sn in c('ssp126','ssp370','ssp585')){
      for(y in c("2021","2041","2061","2081")){
        dat4pre_future <- read.table(paste("/data/fengkai/grassland/UNOISE/future/Clusters_future_predicted",fm,sn,y,standard, "abundance",reg.method,"value.txt",sep="_"),header=T,row.names = 1,sep="\t")
        colnames(dat4pre_future) <- stringr::str_replace_all(colnames(dat4pre_future),"[.]"," ")
        dat4pre_future <- dat4pre_future[,c("long","lat","pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI",
                                            "Low pH","High pH and TN with low SOC","Low SOC","High SOC","High PWetM",
                                            "Low pH high SOC","High MaxTWarmM","Low pH low MAT","High pH","High NDVI")]
        pred.res.cv <- dat4pre_future[which(rowSums(is.na(dat4pre_future))==0),]
        
        for(kcv in 1:cvnum){
          testnum <- neworder[which(rep(1:cvnum,length=nrow(pre.dat))==kcv)]
          traindat <- pre.dat[-testnum,]
          validdat <- pre.dat[testnum,]
          # construct models with training data set
          if(reg.method=="cubist"){
            library(Cubist)
            cub.mod <- cubist(x=traindat[,-ncol(traindat)],y=traindat[,ncol(traindat)],committees = 100)
            cub.mod.cv <- predict(cub.mod, validdat[,-ncol(validdat)], neighbors = 9)
            rmse <- sqrt(mean((cub.mod.cv - validdat[,ncol(validdat)])^2)) ## Test set RMSE
            R2 <- cor(cub.mod.cv, validdat[,ncol(validdat)])^2 ## Test set R^2
            
            sink(paste("/data/fengkai/grassland/UNOISE/future/pH/cv",reg.method,"predicted_current",kcv,cvnum,standard, "abundance_pH_summary.txt",sep="_"))
            print(summary(cub.mod))
            #print(caret::varImp(cub.mod))
            print(paste("RMSE for cross validation: ",rmse,".",sep=""))
            print(paste("R2 for cross validation: ",R2,".",sep=""))
            sink()
            pred.res.cv <- cbind(pred.res.cv,
                                 predict(cub.mod,dat4pre_future[,colnames(traindat[,-ncol(traindat)])],neighbors = 9))
          }else if(reg.method == "randomforest"){
            
          }

        }
        # take average values for each point with 10-cross-validation results 
        # use the standard deviation to indicate the uncertainty
        pred.res.cv$predpHmean <- rowMeans(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)])
        pred.res.cv$predpHsd <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],sd,MARGIN = 1)
        pred.res.cv$predpHcv <- (pred.res.cv$predpHsd/pred.res.cv$predpHmean)*100
        
        pred.res.cv$predpHhighest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],max,MARGIN = 1)
        pred.res.cv$predpHlowest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],min,MARGIN = 1)
        pred.res.cv$predpHRelativeUncertainty_H_L_Range <- (pred.res.cv$predpHhighest-pred.res.cv$predpHlowest)/pred.res.cv$predpHmean*100
        
        pred.res.cv$`predpHu97.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.975)
        pred.res.cv$`predpHl2.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.025)
        pred.res.cv$predpHUncertainty_975_25_Range <- (pred.res.cv$`predpHu97.5`-pred.res.cv$`predpHl2.5`)/pred.res.cv$predpHmean*100
        
        pred.res.cv$predpHdiffmean <-  pred.res.cv$predpHmean - pred.res.cv$pH 
        pred.res.cv$predpHdiffmeanRelaChange <-  pred.res.cv$predpHdiffmean/pred.res.cv$pH *100
        pred.res.cv$predpHdiffsd <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)]-pred.res.cv$pH,sd,MARGIN = 1)
        pred.res.cv$predpHdiffcv <- (pred.res.cv$predpHdiffsd/pred.res.cv$predpHdiffmean)*100
        pred.res.cv$currentpHstate <- sapply(pred.res.cv$pH,function(x){if(x>7){return(1)}else if(x<7){return(-2)}else{return(0)}})
        pred.res.cv$predpHstate <- sapply(pred.res.cv$predpHmean,function(x){if(x>7){return(5)}else if(x<7){return(-4)}else{return(0)}})
        pred.res.cv$predpHchangestate <- pred.res.cv$currentpHstate + pred.res.cv$predpHstate
        #pred.res.cv$predpHchangeAlkineAndAcid <- sapply(pred.res.cv$predpHchangestate,function(x){if(x== -3){return(1)}else if(x==3){return(-1)}else{return(0)}})
        write.table(pred.res.cv,paste("/data/fengkai/grassland/UNOISE/future/pH/future_predicted_pH",fm,sn,y,standard, "abundance",reg.method,"value.txt",sep="_"),sep="\t",col.names=NA,row.names = T,quote=F)
      }
    }
  }
}

# local universal kriging with future scenarios
library(gstat)
library(maptools)
library(raster)
library(stringr)
GC30_grass <- read.table("/data/fengkai/gis/GlobalLand30/GC30_points.txt",header=T,sep="\t")
wd.shape <- readShapePoly("/data/fengkai/gis/WorldMap/ne_10m_admin_0_countries.shp")
dir.sl <- "/data/fengkai/grassland/UNOISE/future/pH"
map.sl <- list.files(dir.sl,"future_predicted_.*_value.txt",include.dirs = T,full.names = TRUE)
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map.tif")
  #newPlot <- str_replace(pred,".txt","_map.png")
  newKriging <- str_replace(pred,".txt","_kriging.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","predpHmean")]
  colnames(data.obs) <- c("X","Y","VALUE")
  library(automap)
  data.obs.auto <- data.obs
  coordinates(data.obs.auto) =~ X+Y
  v <- autofitVariogram(VALUE~1,input_data = data.obs.auto)
  library(parallel)
  nclus = detectCores()-10
  clus <- c(rep("localhost", nclus))
  cl <- makeCluster(clus, type = "SOCK")
  ncell <- nrow(GC30_grass)
  splt = rep(1:nclus, length.out = ncell)
  newdlst = lapply(as.list(1:nclus), function(w) GC30_grass[splt == w,])
  out.clus <- do.call("rbind", parLapply(cl, newdlst,vgmf=v$var_model, olddata=data.obs,function(lst,olddata,vgmf) 
    gstat::krige(formula=VALUE~1,locations=~X+Y,model=vgmf,data=olddata,newdata=lst, nmax=50, nmin=20)
  ))
  stopCluster(cl)
  library(raster)
  predRaster <- rasterFromXYZ(out.clus[,c(1,2,3)])
  projection(predRaster) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(predRaster,newRas,overwrite=TRUE)
  #write.table(out.clus,newKriging,sep="\t",col.names=NA,row.names = T,quote=F)
}

# Uncertainty CV values for pH at future scenario
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map_CV.tif")
  newKriging <- str_replace(pred,".txt","_kriging_CV.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","predpHcv")]
  colnames(data.obs) <- c("X","Y","VALUE")
  if(length(which(is.na(data.obs$VALUE)))>0){
    data.obs<- data.obs[-which(is.na(data.obs$VALUE)),]
  }
  library(automap)
  data.obs.auto <- data.obs
  coordinates(data.obs.auto) =~ X+Y
  v <- autofitVariogram(VALUE~1,input_data = data.obs.auto)
  library(parallel)
  nclus = detectCores()-10
  clus <- c(rep("localhost", nclus))
  cl <- makeCluster(clus, type = "SOCK")
  ncell <- nrow(GC30_grass)
  splt = rep(1:nclus, length.out = ncell)
  newdlst = lapply(as.list(1:nclus), function(w) GC30_grass[splt == w,])
  out.clus <- do.call("rbind", parLapply(cl, newdlst,vgmf=v$var_model, olddata=data.obs,function(lst,olddata,vgmf) 
    gstat::krige(formula=VALUE~1,locations=~X+Y,model=vgmf,data=olddata,newdata=lst, nmax=50, nmin=20)
  ))
  stopCluster(cl)
  library(raster)
  predRaster <- rasterFromXYZ(out.clus[,c(1,2,3)])
  projection(predRaster) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(predRaster,newRas,overwrite=TRUE)
  #write.table(out.clus,newKriging,sep="\t",col.names=NA,row.names = T,quote=F)
}


# make difference between future and current scenarios for mean pH values in raster format
map.ph <- "/data/fengkai/grassland/UNOISE/future/pH/GC30_pH_current_extracted.tif"
dir.future.ph <- "/data/fengkai/grassland/UNOISE/future/pH"
map.fm.ph <- list.files(dir.future.ph,"future_predicted_pH_.*_map.tif",include.dirs = T,full.names = TRUE)
library(doParallel)
cl <- makeCluster(24)
registerDoParallel(cl)
foreach(future=map.fm.ph) %dopar% {
  current=map.ph
  library(stringr)
  newRasMean <- str_replace(future,".tif","_diff.tif")
  library(raster)
  rasFuture <- raster(future)
  rasCurrent <- raster(current)
  system.time(  rasDiff <- overlay(rasFuture,rasCurrent,fun=function(x,y){return(x-y)}) )
  projection(rasDiff) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(rasDiff,newRasMean,overwrite=TRUE)
}
stopCluster(cl)

