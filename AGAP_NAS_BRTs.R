#R code to develop species distribution models for fluvial fishes based on their native ranges using Boosted Regression Trees (BRTs)
#Developed in support of the USGS Aquatic GAP project
#Code developed by Hao Yu & Arthur Cooper, Research Associates, Department of Fisheries and Wildlife, Michigan State University

#Set working directory and remove all objects from the workspace
setwd('../working_directory')
rm(list = ls(all = TRUE))

#Load and attach necessary libraries
library(dismo) #'gbm.step' function to generate BRT models
library(labdsv) #'matrify' function to flip species data table orientation

####################################
#Import fish data and predictor data
####################################
fish<-read.csv("Fish_data.csv", header = T) #input fish data table
fish_name_itis<-fish[,c(6:8)] #subset table that includes fish species ITIS (Integrated Taxonomic Information System) code, common name, and scientific name

#Separate unrestricted data (no sharing restriction/used in BRT model development) from restricted data (cannot be publicly shared/not used in BRT model development)
fish_unrestricted<-fish[is.na(fish$restricted),] #unrestricted fish data
fish_restricted<-fish[!is.na(fish$restricted),] #restricted fish data
#Input predictor variables table for fluvial stream reaches 
predictors_fluvial<-read.csv("Predictors.csv",header=T) 

#######################################
#Import USGS NAS native HUC8 range data
#######################################

#HUC8 range data obtained from USGS Non-indigenous Aquatic Species (NAS) program
HUC8<-read.csv("HUC8_ranges.csv",header=T) #input species HUC8 range data table
HUC8<-HUC8[,-4] #remove unneeded field
#Restrict HUC8 range data to 'native' HUC8s only (HUC8s with 'introduced' status representing non-native range are removed)
HUC8_native<-HUC8[HUC8$OriginStatus=="Native",]

##########################################################
#Import HUC8s for NHDPlusV2 and Fish List
#Combine HUC locations and predictors
##########################################################

#Spatial location table containing HUC8 for each NHDPlusV2 COMID
spatial<-read.csv("HUC8_streams.csv",header=T)

#Pull HUC8 field and merge with predictor table
spatial_HUC8<-spatial[,c(1,16)]
predictors_fluvial_HUC8<-merge(predictors_fluvial,spatial_HUC8,by.x="comid",by.y="COMID")

#Create fish list and associated HUC8 native ranges
fish_list<-read.csv("Fish_list.csv",header=T)
HUC8_native<-HUC8_native[HUC8_native$Scientific_name %in% fish_list$Scientific_name,]

HUC8_native<-HUC8_native[,1:2]

fish_name<-unique(fish_name_itis[fish_name_itis$scientific_name %in% HUC8_native$Scientific_name,])

############################################
#Develop Boosted Regression Tree (BRT) model
############################################

#Remove NAs if present
SpeciesData<-na.omit(subset(fish_unrestricted,select=c(comid_v2,itis_species,sp_count)))

#Changes input species data table orientation from 'stacked' (species records in rows) to species data in individual columns
Species_Matrix<-matrify(SpeciesData)

#Create copy of species matrix that will later be converted from abundances to binary (0/1) presence/absence
Species01<-Species_Matrix
nfhp_Species01<-cbind(rownames(Species01),Species01)
names(nfhp_Species01)[1]<-c("comid")

#Create list of species and native HUC8s
fish_listrange<-merge(fish_name,HUC8_native, by.x="scientific_name",by.y="Scientific_name",all=FALSE)
fish_listrange<-fish_listrange[order(fish_listrange$itis_species),]

#Begin BRT loop
for(i in unique(fish_listrange$itis_species)){
  # Extract scientific and common name
  currentfish<-as.matrix(unique(subset(fish_listrange,itis_species==i,select=c(scientific_name,common_name))))
  scientific<-currentfish[1,1]
  scientific<-gsub(" ", "_", scientific)
  common<-currentfish[1,2]
  common<-gsub(" ", "_", common)

#Pull HUC8 range data for target species
HUC_rangeNATIVE<-subset(fish_listrange,itis_species==i,select=c(common_name,HUC8))
  
#Use HUC8 range to create range-wide predictor variable table 
HUC_predictorsNATIVE<-merge(HUC_rangeNATIVE,predictors_fluvial_HUC8,by="HUC8",all=FALSE)
fish_species01<-subset(nfhp_Species01,select=c("comid",i))

#Use HUC8 range to create predictor variable table for presence-absence locations
SpeciesVariableNATIVE<-merge(fish_species01,HUC_predictorsNATIVE,by="comid",all=FALSE)
  names(SpeciesVariableNATIVE)[2]<-c("Abundance")
SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance>0]<-1
SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance<1]<-0

stat.species<-c()

print(unique(SpeciesVariableNATIVE[4]))
i_variables<-setdiff(c(1:ncol(SpeciesVariableNATIVE)),c(1:4,27))

#Set a seed number so that the process can be repeatable
set.seed(10)
#Check the number of species presences to determine the starting lr ("learning") rate in initial BRT model
np<-length(which(SpeciesVariableNATIVE[,"PA"]==1))#number of presences
nab<-length(which(SpeciesVariableNATIVE[,"PA"]==0))#number of absences
percentp<-np/nrow(SpeciesVariableNATIVE)
if (np < 100){lr<-0.01}else{lr<-0.05}

print(paste("no=",np,",percentp=",percentp,",lr=",lr,sep=" "))
BR<-NULL

#Build BRT model
#Reduce learning rate by half if the best tree model has > 1,000 trees
#Maximum number of trees is set at 10,000
count<-0
while(is.null(BR)){
BR<-gbm.step(data=SpeciesVariableNATIVE, gbm.x=i_variables, gbm.y="PA",
             family = "bernoulli", tree.complexity = 5, learning.rate = lr, max.trees = 10000,
             plot.main=FALSE, keep.fold.models=TRUE, keep.fold.vector=TRUE, keep.fold.fit=TRUE)
       
      
if(!is.null(BR)){
BR_stat<-as.data.frame(cbind(BR$gbm.call$tree.complexity,BR$gbm.call$learning.rate,BR$gbm.call$best.trees))
        names(BR_stat)<-c("tree.complexity","learning.rate","best.n.trees")
print(BR_stat)
if(BR$gbm.call$best.trees<1000) {
BR<-NULL
}
}
      
count<-count+1 # to avoid endless if attempting with too few species presences
if(count>=10){
BR$gbm.call$tree.complexity<-9999
BR$gbm.call$learning.rate<-lr
BR$gbm.call$best.trees<-np
}
lr<-lr*0.5
}
    
dev_exp<-1-(BR$cv.statistics$deviance.mean/BR$self.statistics$mean.null)#model deviance explained

#Gather BRT variable contributions output
varcont<-as.data.frame(BR$contributions)
varcont<-varcont[order(varcont$var),]
varcont<-as.data.frame(t(varcont))
varcont<-varcont[-1,]
varcont["Name"]<-scientific
varcont<-varcont[c(23,1:22)]
    
#Write results of variable contributions to a .csv table
if(i==min(unique(fish_listrange$itis_species))){
write.table(varcont,paste("results/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE)
}else{
write.table(varcont,paste("results/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE,
            col.names=FALSE,append=TRUE)
}

#Compile BRT model statistics
BR_stat<-as.data.frame(cbind(np,nab,percentp,BR$gbm.call$tree.complexity,BR$gbm.call$learning.rate,BR$gbm.call$best.trees,dev_exp,BR$self.statistics$discrimination,BR$cv.statistics$discrimination.mean))
names(BR_stat)<-c("presences", "absences", "prevalence", "tree.complexity","learning.rate","best.n.trees","dev_exp","train.auc","cv.auc")
stat<-cbind(BR_stat)
row.names(stat)<-scientific
stat.species<-rbind(stat.species,stat)

#Write BRT statistics to a new table
if(i==min(unique(fish_listrange$itis_species))){
write.table(stat.species,paste("results/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(stat.species,paste("results/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Develop partial dependence plot of the top 12 predictors
gbm.plot(BR, n.plots=12, write.title = FALSE, common.scale = FALSE,plot.layout=c(3, 4))
  mtext(paste(gsub("_", " ", scientific)," (",gsub("_", " ", common),")",sep=""), outer = TRUE, line=-2, cex = 1.5)
  
savePlot(filename=paste("results/BRT_",scientific,"_",common,"_plots.pdf",sep=""),type=c("pdf"), device=dev.cur())

##############################
#Cross validation of BRT model
##############################    

pred_sum_BRT<-c()
predict_all<-c()
n.fold<-10
k<-0
for(k in 1:n.fold){
selector<-BR$fold.vector
i_fold_t<-which(selector!=k)
i_fold_v<-which(selector==k)
k.cv.tdata<-SpeciesVariableNATIVE[i_fold_t,]
k.cv.vdata<-SpeciesVariableNATIVE[i_fold_v,]

k.cv.fit<-BR$fold.models[[k]]

r_TEST<-predict(k.cv.fit, newdata=k.cv.vdata,n.trees=k.cv.fit$n.trees, type="response")
d_TEST <- as.data.frame(cbind(k.cv.vdata$PA, r_TEST))
dd_TEST<-as.data.frame(cbind(k.cv.vdata, r_TEST))
pred_sum_BRT<-rbind(pred_sum_BRT,d_TEST)
predict_all<-rbind(predict_all,dd_TEST)
}

Deviance_TEST_BRT<-calc.deviance(obs=pred_sum_BRT[,1], pred=pred_sum_BRT[,2], family="bernoulli",calc.mean=TRUE)
pres_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==1, 2]
abs_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==0, 2]

e_TEST_BRT <- evaluate(p=pres_TEST_BRT, a=abs_TEST_BRT)

#Plot and save AUC results
plot(e_TEST_BRT, 'ROC')

savePlot(filename=paste("results/BRT_",scientific,"_",common,"_AUC.tiff",sep=""),type=c("tiff"), device=dev.cur())

########################################################
#Develop presence/absence cutoffs and diagnostic metrics
########################################################

t_TEST_BRT <- threshold(e_TEST_BRT)
cutoff_TEST_BRT<-np/(np+nab)#the prevalence from "threshold" is modeled prevalence

colnames(pred_sum_BRT)<-c("PA","Predict")
pred_sum_BRT<-as.data.frame(pred_sum_BRT)

pred_sum_BRT$predict_PA<-ifelse(pred_sum_BRT$Predict>cutoff_TEST_BRT,1,0)
predict_all$predict_PA<-ifelse(predict_all$r_TEST>cutoff_TEST_BRT,1,0)
names(predict_all)[28]<-c("predict_prob")

conf_d_TEST_BRT<-table(pred_sum_BRT$predict_PA,pred_sum_BRT$PA)

sensitivity_TEST_BRT<-conf_d_TEST_BRT[2,2]/(conf_d_TEST_BRT[2,2]+conf_d_TEST_BRT[2,1])
specificity_TEST_BRT<-conf_d_TEST_BRT[1,1]/(conf_d_TEST_BRT[1,1]+conf_d_TEST_BRT[1,2])

TSS_TEST_BRT<-sensitivity_TEST_BRT+specificity_TEST_BRT-1

det_cv_fold_BRT<-as.data.frame(cbind(Deviance_TEST_BRT, np, nab,e_TEST_BRT@auc,e_TEST_BRT@cor,cutoff_TEST_BRT,sensitivity_TEST_BRT,specificity_TEST_BRT,TSS_TEST_BRT))

colnames(det_cv_fold_BRT)<-c("Deviance","np","na","auc","cor","threshold","sensitivity","specificity","TSS")
rownames(det_cv_fold_BRT)<-scientific
if(i==min(unique(fish_listrange$itis_species))){
write.table(det_cv_fold_BRT,paste("results/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(det_cv_fold_BRT,paste("results/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Write output cross validation tables
write.csv(pred_sum_BRT,paste("results/BRT_CV_predict_",scientific,"_",common,".csv"))
write.csv(predict_all,paste("results/BRT_CV_predict_all_",scientific,"_",common,".csv"))


########################################################################
#Project model results to all fluvial stream reaches within native range
########################################################################

predictor_species<-HUC_predictorsNATIVE[,c(4:25)]
predict_native_region<-predict(BR,predictor_species,n.trees=BR$gbm.call$best.trees,type="response" )

predict_native_region_prob<-as.data.frame(cbind(HUC_predictorsNATIVE[,c(1:3)],predict_native_region))
names(predict_native_region_prob)[4]<-c("predict_prob")

predict_native_region_prob$predict_PA<-ifelse(predict_native_region_prob$predict_prob>cutoff_TEST_BRT,1,0)

write.csv(predict_native_region_prob,paste("results/BRT_",scientific,"_",common,"_native_prediction.csv",sep=""))

save(BR,file=paste("results/BRT_",scientific,"_",common,"_.RData",sep=""))# Save to R data file
}#BRT model loop

#END 
################

