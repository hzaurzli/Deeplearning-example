
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

env_GC30_grass_current <- read.table("/data/fengkai/gis/GlobalLand30/GC30Env_current_extracted.txt",header=T,row.names = 1,sep="\t")
if(!dir.exists("/data/fengkai/grassland/UNOISE/current")){dir.create("/data/fengkai/grassland/UNOISE/current")}
standard = "Relative"
cvnum = 10 #　10-fold cross validation
for(reg.method in c("cubist")){ # "randomforest",
  for (slc in eco) {
    otu.sub <- OTUcoreS1[eco.grp[[slc]],]
    pre.dat <- cbind(env.sd[,c("pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")],colSums(otu.sub))
    if(standard == "Standardized"){
      pre.dat[,ncol(pre.dat)] <- decostand(pre.dat[,ncol(pre.dat)],method = "standardize") 
    }
    colnames(pre.dat)[ncol(pre.dat)] <- "abundance"
    set.seed(36)
    neworder <- sample(1:nrow(pre.dat),nrow(pre.dat),replace = FALSE)
    dat4pre_now <- env_GC30_grass_current[,c("long","lat","soil_phh2o_0_15","soil_soc_0_15","soil_tn_0_15","BIO1","BIO5","BIO12","BIO13","plant_ndvi_2010_2020")]
    colnames(dat4pre_now) <- c("long","lat","pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")
    pred.res.cv <- dat4pre_now[which(rowSums(is.na(dat4pre_now))==0),]
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
        
        sink(paste("/data/fengkai/grassland/UNOISE/current/cv",reg.method,"predicted_current",kcv,cvnum,standard, "abundance",slc,"summary.txt",sep="_"))
        print(summary(cub.mod))
        print(caret::varImp(cub.mod))
        print(paste("RMSE for cross validation: ",rmse,".",sep=""))
        print(paste("R2 for cross validation: ",R2,".",sep=""))
        sink()
        pred.res.cv <- cbind(pred.res.cv,
                         predict(cub.mod,dat4pre_now[which(rowSums(is.na(dat4pre_now))==0),-c(1,2)],neighbors = 9))
      }else if(reg.method == "randomforest"){
        
      }
    }
    # take average values for each point with 10-cross-validation results 
    # use the standard deviation to indicate the uncertainty
    pred.res.cv$mean <- rowMeans(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)])
    pred.res.cv$sd <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],sd,MARGIN = 1)
    pred.res.cv$cv <- (pred.res.cv$sd/pred.res.cv$mean)*100
    
    pred.res.cv$highest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],max,MARGIN = 1)
    pred.res.cv$lowest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],min,MARGIN = 1)
    pred.res.cv$RelativeUncertainty_H_L_Range <- (pred.res.cv$highest-pred.res.cv$lowest)/pred.res.cv$mean*100
    
    pred.res.cv$`u97.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.975)
    pred.res.cv$`l2.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.025)
    pred.res.cv$Uncertainty_975_25_Range <- (pred.res.cv$`u97.5`-pred.res.cv$`l2.5`)/pred.res.cv$mean*100
    write.table(pred.res.cv,paste("/data/fengkai/grassland/UNOISE/current/current_predicted",standard, "abundance",slc,reg.method,"value.txt",sep="_"),sep="\t",col.names=NA,row.names = T,quote=F)
  }
}

# construct cubist model and predict the values to the future map
env_GC30_grass_current <- read.table("/data/fengkai/gis/GlobalLand30/GC30Env_current_extracted.txt",header=T,row.names = 1,sep="\t")
if(!dir.exists("/data/fengkai/grassland/UNOISE/future")){dir.create("/data/fengkai/grassland/UNOISE/future")}
standard = "Relative"
cvnum = 10 #　10-fold cross validation
library(doParallel)
cl <- makeCluster(11)
registerDoParallel(cl)
for(fm in c("BCC_CSM2_MR","IPSL_CM6A_LR")){
  for(sn in c('ssp126','ssp370','ssp585')){
    for(y in c("2021","2041","2061","2081")){
      biofu <- read.table(paste("/data/fengkai/gis/GlobalLand30/GC30Env",fm,sn,y,"BIO1_19.txt",sep="_"),header=T,row.names = 1,sep="\t")
      for(reg.method in c("cubist")){ # "randomforest",
        foreach (slc=eco) %dopar%  {
          otu.sub <- OTUcoreS1[eco.grp[[slc]],]
          pre.dat <- cbind(env.sd[,c("pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")],colSums(otu.sub))
          if(standard == "Standardized"){
            pre.dat[,ncol(pre.dat)] <- decostand(pre.dat[,ncol(pre.dat)],method = "standardize") 
          }
          colnames(pre.dat)[ncol(pre.dat)] <- "abundance"
          set.seed(36)
          neworder <- sample(1:nrow(pre.dat),nrow(pre.dat),replace = FALSE)
          dat4pre_future <- cbind(env_GC30_grass_current[,c("long","lat","soil_phh2o_0_15","soil_soc_0_15","soil_tn_0_15")],
                                  biofu[,c("FBIO1","FBIO5","FBIO12","FBIO13")],env_GC30_grass_current[,c("plant_ndvi_2010_2020")])
          colnames(dat4pre_future) <- c("long","lat","pH","SOC","TN","MAT","MaxTWarmM","MAP","PWetM","NDVI")
          pred.res.cv <- dat4pre_future[which(rowSums(is.na(dat4pre_future))==0),]
          for(kcv in 1:cvnum){
            testnum <- neworder[which(rep(1:cvnum,length=nrow(pre.dat))==kcv)]
            traindat <- pre.dat[-testnum,]
            validdat <- pre.dat[testnum,]
            # construct models with training data set
            if(reg.method=="cubist"){
              library(Cubist)
              cub.mod <- cubist(x=traindat[,-ncol(traindat)],y=traindat[,ncol(traindat)],committees = 100)
              pred.res.cv <- cbind(pred.res.cv,
                                   predict(cub.mod,dat4pre_future[which(rowSums(is.na(dat4pre_future))==0),-c(1,2)],neighbors = 9))
            }else if(reg.method == "randomforest"){
              
            }
          }
          # take average values for each point with 10-cross-validation results 
          # use the standard deviation to indicate the uncertainty
          pred.res.cv$mean <- rowMeans(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)])
          pred.res.cv$sd <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],sd,MARGIN = 1)
          pred.res.cv$cv <- (pred.res.cv$sd/pred.res.cv$mean)*100
          
          pred.res.cv$highest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],max,MARGIN = 1)
          pred.res.cv$lowest <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],min,MARGIN = 1)
          pred.res.cv$RelativeUncertainty_H_L_Range <- (pred.res.cv$highest-pred.res.cv$lowest)/pred.res.cv$mean*100
          
          pred.res.cv$`u97.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.975)
          pred.res.cv$`l2.5` <- apply(pred.res.cv[,grep("predict",colnames(pred.res.cv),perl=T)],quantile,MARGIN = 1,probs=0.025)
          pred.res.cv$Uncertainty_975_25_Range <- (pred.res.cv$`u97.5`-pred.res.cv$`l2.5`)/pred.res.cv$mean*100
          write.table(pred.res.cv,paste("/data/fengkai/grassland/UNOISE/future/future_predicted",fm,sn,y,standard, "abundance",slc,reg.method,"value.txt",sep="_"),sep="\t",col.names=NA,row.names = T,quote=F)
        }
      }
    }
  }
}
stopCluster(cl)

# local universal kriging with current scenarios
library(gstat)
library(maptools)
library(raster)
GC30_grass <- read.table("/data/fengkai/gis/GlobalLand30/GC30_points.txt",header=T,sep="\t")
wd.shape <- readShapePoly("/data/fengkai/gis/WorldMap/ne_10m_admin_0_countries.shp")
dir.sl <- "/data/fengkai/grassland/UNOISE/current"
map.sl <- list.files(dir.sl,"current_predicted_.*_value.txt",include.dirs = T,full.names = TRUE)
library(stringr)
# mean values for each cluster at current scenario
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map.tif")
  #newPlot <- str_replace(pred,".txt","_map.png")
  newKriging <- str_replace(pred,".txt","_kriging.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","mean")]
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

# Uncertainty CV values for each cluster at current scenario
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map_CV.tif")
  newKriging <- str_replace(pred,".txt","_kriging_CV.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","cv")]
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

# 
# local universal kriging with future scenarios
library(gstat)
library(maptools)
library(raster)
library(stringr)
GC30_grass <- read.table("/data/fengkai/gis/GlobalLand30/GC30_points.txt",header=T,sep="\t")
wd.shape <- readShapePoly("/data/fengkai/gis/WorldMap/ne_10m_admin_0_countries.shp")
dir.sl <- "/data/fengkai/grassland/UNOISE/future"
map.sl <- list.files(dir.sl,"future_predicted_.*_value.txt",include.dirs = T,full.names = TRUE)
# map.sl <- list.files(dir.sl,"future_predicted_BCC_CSM2_MR_.*_value.txt",include.dirs = T,full.names = TRUE)
# [c(14,16,18,19,20,34,36,38,39,40,94,96,98,99,100,114,116,118,119,120)]
# [(1:120)[-c(14,16,18,19,20,34,36,38,39,40,94,96,98,99,100,114,116,118,119,120)]]
# mean values for each cluster at future scenario
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map.tif")
  #newPlot <- str_replace(pred,".txt","_map.png")
  newKriging <- str_replace(pred,".txt","_kriging.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","mean")]
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

# Uncertainty CV values for each cluster at future scenario
for(pred in map.sl){
  newRas <- str_replace(pred,".txt","_map_CV.tif")
  newKriging <- str_replace(pred,".txt","_kriging_CV.txt")
  pred.res.cv <- read.table(pred,header=T,row.names = 1,sep="\t")
  data.obs <- pred.res.cv[,c("long","lat","cv")]
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

# make difference between future and current scenarios for mean and CV values in raster format
dir.c <- "/data/fengkai/grassland/UNOISE/current"
map.cm <- list.files(dir.c,"current_predicted_.*_map.tif",include.dirs = T,full.names = TRUE)
#map.ccv <- list.files(dir.c,"current_predicted_.*_map_CV.tif",include.dirs = T,full.names = TRUE)
dir.f <- "/data/fengkai/grassland/UNOISE/future"
map.fm <- list.files(dir.f,"future_predicted_BCC_CSM2_.*_map.tif",include.dirs = T,full.names = TRUE)
#map.fm <- list.files(dir.f,"future_predicted_IPSL_CM6A_LR_.*_map.tif",include.dirs = T,full.names = TRUE)
#map.fcv <- list.files(dir.f,"future_predicted_BCC_CSM2_.*_map_CV.tif",include.dirs = T,full.names = TRUE)
library(doParallel)
cl <- makeCluster(40)
registerDoParallel(cl)
foreach(future=map.fm) %dopar% {
  library(stringr)
  newRasMean <- str_replace(future,".tif","_diff.tif")
  newRasCV <- str_replace(future,".tif","_diff_CV.tif")
  future.cv <- str_replace(future,".tif","_CV.tif")
  ecoCluster <- unlist(str_split(future,"_"))[10]
  current <- map.cm[grep(paste("",ecoCluster,"",sep="_"),map.cm)]
  current.cv <- str_replace(current,".tif","_CV.tif")
  library(raster)
  rasFuture <- raster(future)
  rasFutureCV <- raster(future.cv)
  rasCurrent <- raster(current)
  rasCurrentCV <- raster(current.cv)
  system.time(  rasDiff <- overlay(rasFuture,rasCurrent,fun=function(x,y){return(x-y)}) )
  projection(rasDiff) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(rasDiff,newRasMean,overwrite=TRUE)
  system.time( rasDiffCV <- overlay(rasFuture,rasFutureCV,rasCurrent,rasCurrentCV,fun=function(v1,cv1,v2,cv2){
    x <- (v1*cv1+v2*cv2)/abs(v1-v2)
    x[is.infinite(x)] <- 0
    x[is.nan(x)] <-0
    x[is.na(x)] <- 0
    return(x)
  }   )  )
  projection(rasDiffCV) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  writeRaster(rasDiffCV,newRasCV,overwrite=TRUE)
}
stopCluster(cl)
