
######################################################
###                                                ### 
###       Calibrate RF algo releves ~ layers       ###
###           projet ORION - avril 2022            ### 
###                                                ###
######################################################



rm(list=ls())


library(raster); library(rgeos); library(rgdal)
library(pdp);library(caret);library(randomForest);library(TeachingDemos)



mask=raster::mask
extract=raster::extract
crop=raster::crop
writeRaster=raster::writeRaster

#rasterOptions(tmpdir="D:/Rtemp/")
options(digits = 15)


#setwd("Y:/ORION/")

load("_data/19042022_RF.RData")


# Write table and rasters -------------------------------------------------

Calib = read.csv("CALIB/EM_variables_CC.csv")[,-1]


DEM_CCVCMB <- raster("CALIB/Variables_explicatives/CCVCMB/DEM_CC.tif") 
TPI_CCVCMB <- raster("CALIB/Variables_explicatives/CCVCMB/TPI_CC.tif") 
DAH_CCVCMB <- raster("CALIB/Variables_explicatives/CCVCMB/DAH_CC.tif") 
SFGDD_median_CC <- raster("CALIB/Variables_explicatives/CCVCMB/SFGDD_CC.tif") 
P_ampli_CC <- raster("CALIB/Variables_explicatives/CCVCMB/Ampli_CC.tif") 
P_EOSD_CC <- raster("CALIB/Variables_explicatives/CCVCMB/EOSD_CC.tif")
P_MAXD_CC <- raster("CALIB/Variables_explicatives/CCVCMB/MAXD_CC.tif")
P_MAXV_CC <- raster("CALIB/Variables_explicatives/CCVCMB/MAXV_CC.tif") 
P_NARI_CC <- raster("CALIB/Variables_explicatives/CCVCMB/NARI_CC.tif") 
P_SOSD_CC <- raster("CALIB/Variables_explicatives/CCVCMB/SOSD_CC.tif") 
P_SPROD_CC <- raster("CALIB/Variables_explicatives/CCVCMB/SPROD_CC.tif") 




# test correlations between variables -------------------------------------

# for phenological variables in non-vegetated contexts, replace NA w/ 0 - double check this w/ Arthur

Calib[which(is.na(Calib$AMPLI) | is.na(Calib$SOSD)), c(9:12,14:15)] <- 0


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(Calib[,c(5:15)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)



# RM AMPLI, EOSD, MAXD & MAXV (too correlated with other variables)

Calib <- Calib[,-c(5, 9,10,11, 12)]


pairs(Calib[,c(5:10)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)





# calibrate RF algorithm (based on Arthur's script) --------------------------------------------------


# create first dataset without under-represented classes (Calib 1)

Calib1 <- Calib[-which(Calib$Class %in% c("Rhodo", "Juniper", "Water", "Alpine_grassland")),]
Calib1 <- Calib1[,-c(1:3)] # rm ID and XY coordinates

Calib1$Class <- as.factor(Calib1$Class)


### fit model
control <- trainControl(method="oob", number=20, search="random")
rf.trained <- train(
  Class~., data=Calib1, method="rf", metric="Accuracy", trControl=control, ntree = 2000, importance = T)

# Variable importance
print('estimating variable importance...')
rf.best.fit <- rf.trained$finalModel
var.imp.dt <- rf.best.fit$importance

# Extract accuracy information
DATA.eval$predicted <- predict.train(rf.trained, newdata = DATA.eval[,-1])
conf.mtx <- confusionMatrix(reference = DATA.eval$CLASS, data = DATA.eval$predicted, mode='sens_spec')
Overall.Accur[i] = mean(conf.mtx$byClass[,11])


# Partial dependencies

#lande
V1_SFGDD_Lande <- partial(object = rf.trained, progress = "text", pred.var = "SFGDD", which.class = "Lande", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_NARI_Lande <- partial(object = rf.trained, progress = "text", pred.var = "NARI", which.class = "Lande", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SOSD_Lande <- partial(object = rf.trained, progress = "text", pred.var = "SOSD", which.class = "Lande", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SPROD_Lande <- partial(object = rf.trained, progress = "text", pred.var = "SPROD", which.class = "Lande", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)


#prairie
V1_SFGDD_Grass <- partial(object = rf.trained, progress = "text", pred.var = "SFGDD", which.class = "Nardus_grassland", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_NARI_Grass <- partial(object = rf.trained, progress = "text", pred.var = "NARI", which.class = "Nardus_grassland", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SOSD_Grass <- partial(object = rf.trained, progress = "text", pred.var = "SOSD", which.class = "Nardus_grassland", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SPROD_Grass <- partial(object = rf.trained, progress = "text", pred.var = "SPROD", which.class = "Nardus_grassland", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)

#tall_shrub
V1_SFGDD_TS <- partial(object = rf.trained, progress = "text", pred.var = "SFGDD", which.class = "Tall_shrub", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_NARI_TS <- partial(object = rf.trained, progress = "text", pred.var = "NARI", which.class = "Tall_shrub", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SOSD_TS <- partial(object = rf.trained, progress = "text", pred.var = "SOSD", which.class = "Tall_shrub", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)
V1_SPROD_TS <- partial(object = rf.trained, progress = "text", pred.var = "SPROD", which.class = "Tall_shrub", train = Calib1, smooth = T, type = 'classification', prob = T, trim.outliers = T)



### Plot response curves

par(mfrow=c(2,2))
plot(V1_SFGDD_Lande$SFGDD, V1_SFGDD_Lande$yhat, type='l',lwd=2, xlab="SFGDD", ylim=range(0,0.4), ylab="Classification probability", col="red")
points(V1_SFGDD_Grass$SFGDD, V1_SFGDD_Grass$yhat, type='l', lwd=2, col="blue")
points(V1_SFGDD_TS$SFGDD, V1_SFGDD_TS$yhat, type='l', lwd=2, col="green")

plot(V1_NARI_Lande$NARI, V1_NARI_Lande$yhat, type='l',lwd=2, ylim=range(0,0.4),xlab="NARI", ylab="Classification probability", col="red")
points(V1_NARI_Grass$NARI, V1_NARI_Grass$yhat, type='l', lwd=2,col="blue")
points(V1_NARI_TS$NARI, V1_NARI_TS$yhat, type='l', lwd=2,col="green")

plot(V1_SOSD_Lande$SOSD, V1_SOSD_Lande$yhat, type='l', lwd=2, ylim=range(0,0.4), xlab="SOSD", ylab="Classification probability", col="red")
points(V1_SOSD_Grass$SOSD, V1_SOSD_Grass$yhat, type='l', lwd=2, col="blue")
points(V1_SOSD_TS$SOSD, V1_SOSD_TS$yhat, type='l', lwd=2, col="green")

plot(V1_SPROD_Lande$SPROD, V1_SPROD_Lande$yhat, type='l', lwd=2, ylim=range(0,0.4), xlab="SPROD", ylab="Classification probability", col="red")
points(V1_SPROD_Grass$SPROD, V1_SPROD_Grass$yhat, type='l', lwd=2,col="blue")
points(V1_SPROD_TS$SPROD, V1_SPROD_TS$yhat, type='l', lwd=2,col="green")


legend(50, 0.4, lwd=2, col=c("red", "blue", "green"), legend=c("Lande", "Prairie", "Tall-shrub"), cex=0.75)





# PREDICT : create new data -----------------------------------------------


# NEWDATA from rasters

P_SOSD_CC[which(is.na(P_SOSD_CC[]))] <- 0
P_SOSD_CC <- mask(P_SOSD_CC, SFGDD_median_CC)

P_SPROD_CC[which(is.na(P_SPROD_CC[]))] <- 0
P_SPROD_CC <- mask(P_SPROD_CC, SFGDD_median_CC)

Var_exp <- data.frame(TPI = TPI_CCVCMB[which(!is.na(TPI_CCVCMB[]))],
                      DAH = DAH_CCVCMB[which(!is.na(DAH_CCVCMB[]))],
                      SFGDD = SFGDD_median_CC[which(!is.na(SFGDD_median_CC[]))],
                      NARI = P_NARI_CC[which(!is.na(P_NARI_CC[]))],
                      SOSD = P_SOSD_CC[which(!is.na(P_SOSD_CC[]))],
                      SPROD = P_SPROD_CC[which(!is.na(P_SPROD_CC[]))])

# PREDICT

Var_exp$predicted <- predict.train(rf.trained, newdata = Var_exp)


# Assign code to predicted classes
Var_exp$Class_pred <- rep(NA, nrow(Var_exp))

# 0 == Urban
# 1 == Montane grassland
# 2 == Forest
# 3 == Tall shrub
# 4 == Lande
# 5 == Prairie (? nard)
# 6 == Rock
# 7 == Glacier


Var_exp[which(Var_exp$predicted=="Urban"), "Class_pred"] <- 0
Var_exp[which(Var_exp$predicted=="Montane_grassland"), "Class_pred"] <- 1
Var_exp[which(Var_exp$predicted=="Forest"), "Class_pred"] <- 2
Var_exp[which(Var_exp$predicted=="Tall_shrub"), "Class_pred"] <- 3
Var_exp[which(Var_exp$predicted=="Lande"), "Class_pred"] <- 4
Var_exp[which(Var_exp$predicted=="Nardus_grassland"), "Class_pred"] <- 5
Var_exp[which(Var_exp$predicted=="Rock"), "Class_pred"] <- 6
Var_exp[which(Var_exp$predicted=="Glacier"), "Class_pred"] <- 7





# store predicted values in raster ----------------------------------------

LC_pred <- DEM_CCVCMB
LC_pred[which(!is.na(LC_pred[]))] <- Var_exp$Class_pred


plot(LC_pred)


writeRaster(LC_pred, "Y:/ORION/PRED/LC_pred_7_19042022.tif")


Forest <- reclassify(LC_pred, c(-Inf, 1.5, NA,  1.5, 2.5, 1,  2.5, Inf, NA))
TS <- reclassify(LC_pred, c(-Inf, 2.5, NA,  2.5, 3.5, 1,  3.5, Inf, NA))
Lande <- reclassify(LC_pred, c(-Inf, 3.5, NA,  3.5, 4.5, 1, 4.5, Inf, NA))
Grass <- reclassify(LC_pred, c(-Inf, 4.5, NA,  4.5, 5.5, 1, 5.5, Inf, NA))


writeRaster(Forest, "Y:/ORION/PRED/Forest_19042022.tif")
writeRaster(TS, "Y:/ORION/PRED/Tallshrub_19042022.tif")
writeRaster(Lande, "Y:/ORION/PRED/Lande_19042022.tif")
writeRaster(Grass, "Y:/ORION/PRED/Grass_19042022.tif")






# # create second dataset with simplified classificaiton (Calib 2)
# 
# Calib2 <- Calib[,-c(1:3)] # rm ID and XY coordinates
# 
# Calib2[which(Calib2$Class == "Rhodo"), "Class"] <- "Lande"
# Calib2[which(Calib2$Class == "Juniper"), "Class"] <- "Lande"
# Calib2[which(Calib2$Class == "Alpine_grassland"), "Class"] <- "Nardus_grassland"
# 
# Calib2$Class <- as.factor(Calib2$Class)
# 
# 
# ### fit model
# control <- trainControl(method="oob", number=20, search="random")
# rf.trained <- train(
#   Class~., data=Calib2, method="rf", metric="Accuracy", trControl=control, ntree = 2000, importance = T)
# 
# # Variable importance
# print('estimating variable importance...')
# rf.best.fit <- rf.trained$finalModel
# var.imp.dt <- rf.best.fit$importance
# 
# # Extract accuracy information
# DATA.eval$predicted <- predict.train(rf.trained, newdata = DATA.eval[,-1])
# conf.mtx <- confusionMatrix(reference = DATA.eval$CLASS, data = DATA.eval$predicted, mode='sens_spec')
# Overall.Accur[i] = mean(conf.mtx$byClass[,11])




























# ARTHUR'S original script ------------------------------------------------


# Load data
DATA = read.table("D:/Recherche/Publications/2021/Glacier Noir/DATA/TABLE/DATA_NDVI_UAV_GN_ALL_DATA.csv",sep=",",header=T)
TRASP = read.table("D:/Recherche/Publications/2021/Glacier Noir/DATA/TABLE/DATA_NDVI_UAV_GN_ALL_DATA_TRASP.csv",sep=",",header=T)

# Prepare classes
DATA$CLASS = NA
DATA$CLASS[which(DATA$RES > 0.025)] = "Pos"
DATA$CLASS[which(DATA$RES < -0.025)] = "Neg"
DATA$CLASS[which(DATA$RES < 0.025 & DATA$RES > -0.025)] = "Nul"
table(DATA$CLASS)

# Equalize samples size
DATA = DATA[-sample(which(DATA$CLASS == "Neg"),size = 15),]
DATA = DATA[-sample(which(DATA$CLASS == "Nul"),size = 9),]

# Select variables to use in RF
DATA = DATA[,c("CLASS","VEG_HEIGHT_MEAN","TRASP_MEAN","TRASP_STD","SLOPE_MEAN","SLOPE_STD")]
DATA$CLASS = as.factor(DATA$CLASS)

# Seek for highly correlated variables to remove (r? > 0.75)
print('seek for highly corrected variables...')
high.cors = cor(DATA[,sapply(DATA,is.numeric)],use = "pair")
if(any(high.cors > 0.75 & high.cors != 1) |
   any(high.cors < -0.75 & high.cors != 1)){print('find one !')} else {print('None...')}

# Dispatch 2/3 between training and evaluation
print('dispatching data 2/3...')
train.indx <- sample(1:nrow(DATA), size = nrow(DATA)*0.666)
DATA.train <- DATA[train.indx,]
DATA.eval <- DATA[-train.indx,]

# Fitting model
print('fitting models...')
control <- trainControl(method="oob", number=20, search="random")
rf.trained <- train(CLASS~., data=DATA.train, method="rf", metric="Accuracy", trControl=control, ntree = 2000, importance = T)

# Variable importance
print('estimating variable importance...')
rf.best.fit <- rf.trained$finalModel
var.imp.dt <- rf.best.fit$importance

# Extract accuracy information
DATA.eval$predicted <- predict.train(rf.trained, newdata = DATA.eval[,-1])
conf.mtx <- confusionMatrix(reference = DATA.eval$CLASS, data = DATA.eval$predicted, mode='sens_spec')
Overall.Accur[i] = mean(conf.mtx$byClass[,11])

# Partial dependencies
V1_pos_pdp <- partial(object = rf.trained, progress = "text", pred.var = "VEG_HEIGHT_MEAN", which.class = "Pos", train = DATA, smooth = T, type = 'classification', prob = T, trim.outliers = T)



