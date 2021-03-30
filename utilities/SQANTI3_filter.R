# Authors: Cecile Pereira, Lorena de la Fuente, Francisco Pardo, Leandro Balzano-Nogueira
# f.pardo.palacios@gmail.com
# Universitat Politècnica de València (UPV)
# leobalzano@ufl.edu
# Genetics Institute, University of Florida (Gainesville)
# Last update: Mar/30/2021

###############################
##### Package requirements ####
###############################

list.of.packages <- c("rpart", "ROCR", "caret", "optparse", "ggplot2", "lattice", "foreach", "e1071","randomForest", "partykit", "ipred", "rpart.plot", "doMC", "nnet", "ROSE", "pROC", "MLmetrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

rm(list=ls())
library(rpart)           #Rpart package
library(partykit)        #"Elegant" tree design package
library(ipred)           #Improved Predictors package
library(rpart.plot)      #Tree design 
library(ROCR)            #ROC curve
library(caret)           #confusion matrix: caret_6.0-76 Version
library(randomForest)    #randomForest
library(nnet)            #neural networks
library("ggplot2")
#library('RWeka')
library('ROSE')
library("optparse")
library(pROC)
library('MLmetrics')

###############################
########## Functions ##########
###############################
#input: result of rfe function
filter_wilcox <-function(profile.1){
  # compare mean values of KAPPA in function of the number of variables:
  wilcox=pairwise.wilcox.test(profile.1$resample$Kappa, profile.1$resample$Variables,p.adjust.method='none') #<0.1: selection of 4 features
  ttest=pairwise.t.test(profile.1$resample$Kappa, profile.1$resample$Variables,p.adjust.method='none') #<0.1: selection of 4 features
  
  #first column with all values > 0.1 is the min number of features.
  minvalues=apply(na.omit(wilcox$p.value),2,min)
  which(minvalues>0.1)
  numbervariable=which(minvalues>0.1)[1]
  
  print("Features selected by the recursive feature elimination")
  print(profile.1$optVariables[1:numbervariable])
  return(numbervariable)
}

###############################
### Definition of the inputs ##
###############################

option_list=list(
  make_option(c("-c","--sqanti_classif"),type="character", default = NULL, help="Sqanti classification output file"),
  make_option(c("-d","--dir"),type="character", default="Filter_out", help="Output directory name"),
  make_option(c('-w','--wilcox'),type="integer",default=0, help="Feature selection option: 0 = max rfe value, 1 = minimum number of feature with a mean non different to the mean of the classifier obtain on all the features."),
  make_option(c("-t","--pourcent_training"),type="double",default=0.8,help="the percentage of data that goes to training (parameter p of the function createDataPartition)"),
  make_option(c("-e","--exon_cov"),type="integer", default=0,help="Should the exon coverage be considered ? 0 No, 1 Yes, default 0"),
  make_option(c("-p","--TP"), type="character",default = NULL,help="file containing the list of the TP transcripts, one ID by line, no header"),
  make_option(c("-n","--TN"), type="character",default = NULL,help="file containing the list of the TN transcripts, one ID by line, no header"),
  make_option(c("-j","--threshold"), default = "0.5", help="machine learning probability threshold to classify posiive isoforms"),
  make_option(c("-i","--intrapriming"), default = "80", help="polyA percentage thereshold to flag an isoform as intra-priming")
)
 
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args
opt$threshold = as.numeric(opt$threshold)
opt$intrapriming = as.numeric(opt$intrapriming)

if(opt$TP=="NULL"){
  opt$TP = NULL
}

if(opt$TN=="NULL"){
  opt$TN = NULL
}
 
d <-read.table(file=opt$sqanti_classif, sep ="\t", header=TRUE,as.is =TRUE)
dme <- d[which(d$exons==1),] # store mono-exon transcript for intrapriming filtering.
d1 <- d[which(d$exons!=1),] #remove monoexon

if(nrow(d1)==0){
  print("Machine Learning filtering won't be applied because all the isoforms are mono-exon")
}else{
  print("Number of of transcripts with more than one exon:")
  print(dim(d1)[1])
}

if(is.null(opt$TP)|is.null(opt$TN)){
  print("Training set not provided, It will be created from input data.")
  if(nrow(d1[d1$structural_category=="full-splice_match",])<40 | nrow(d1[which(d1$structural_category=="novel_not_in_catalog" & d1$all_canonical=="non_canonical"),])<40){
        print("Not enought Full splice match (FSM) and/or novel not in catalog-noncanonical (NNC-NC) isoforms to be used as training set. Skipping machine learning.")
        flag=FALSE
  }else{print("Full splice match as true positives, novel not in catalogue as true negatives")
        flag=TRUE}
}

d1$FSM_class=as.factor(d1$FSM_class)
d1$coding=as.factor(d1$coding)
d1$bite=as.factor(d1$bite)
d1$within_cage_peak=as.factor(d1$within_cage_peak)

if (any(grep('^FL\\.', names(d1)))){d1$FL= rowSums(d1[,grep('^FL\\.', names(d1))])}

dtotal=d1

dme$FSM_class=as.factor(dme$FSM_class)
dme$coding=as.factor(dme$coding)
dme$bite=as.factor(dme$bite)
dme$within_cage_peak=as.factor(dme$within_cage_peak)

if (any(grep('^FL\\.', names(dme)))){dme$FL= rowSums(dme[,grep('^FL\\.', names(dme))])}

###############################
######## ML CLASSIFIER ########
###############################

if(flag==TRUE){
  
  print("Data description:")
  
  #############################
  # Remove descriptive columns
  #############################
  
  if(is.null(d1$isoform)){
    rownames(d1)=d1$PacBio
    d1=d1[,-which(colnames(d1)=='PacBio')]
  }else{
    rownames(d1)=d1$isoform
    d1=d1[,-which(colnames(d1)=='isoform')]
  }
  
  ###  línea de Fran
  FL_samples_cols=colnames(d1)[grep('^FL\\.', names(d1))]
  colRem = c('chrom','strand','associated_gene', 'associated_transcript', 'ref_length','ref_exons','diff_to_TSS',
             'diff_to_TTS', 'ORF_length', 'CDS_length', 'CDS_start', 'CDS_end', 'CDS_genomic_start',
             'CDS_genomic_end', 'polyA_dist', 'min_cov_pos', 'seq_A_downstream_TTS')
  d3=d1[,-which(colnames(d1)%in%colRem)]
  #colnames(d3)
  
  
  
  ###############################
  ##### 1) Dealing with NA's ####
  ###############################
  
  ### OJO ALGUNOS TRANSCRITOS EN EL CLASSIFICATION DE SQANTI2 PUEDEN TENER NA DONDE NO DEBERÍA, HAY QUE MIRAR EN EL SCRIPT DE SQANTI2 QC 
  ### EJEMPLO NA EN LA COLUMNA ALL_CANONICAL (EN TRANSCRITOS ANTISENSE). Son los mismos que tienen NA en el bite
  ### polyA_motif también
  ### polyA_dist
  ### dist_to_cage_peak
  ### predicted_NMD
  
  # if NA in polyA motif, replace it by "No_polyA_motif" and if not change it by "Motif_found"
  d3[which(!is.na(d3$polyA_motif)), "found_polyA"] <- "Motif_found"
  d3[which(is.na(d3$polyA_motif)), "found_polyA"] <- "Not_found"
  d3$found_polyA=as.factor(d3$found_polyA)
  d3=d3[,-which(colnames(d3)=="polyA_motif")]
  
  dtotal[which(!is.na(dtotal$polyA_motif)), "found_polyA"] <- "Motif_found"
  dtotal[which(is.na(dtotal$polyA_motif)), "found_polyA"] <- "Not_found"
  dtotal$found_polyA=as.factor(dtotal$found_polyA)
  
  dme[which(!is.na(dme$polyA_motif)), "found_polyA"] <- "Motif_found"
  dme[which(is.na(dme$polyA_motif)), "found_polyA"] <- "Not_found"
  dme$found_polyA=as.factor(dme$found_polyA)
  
  d3$within_assoc_TSS <- FALSE
  d3[which(abs(as.numeric(d3$diff_to_gene_TSS))<50), "within_assoc_TSS"] <- TRUE
  d3$within_assoc_TSS=as.factor(d3$within_assoc_TSS)
  d3=d3[,-which(colnames(d3)=="diff_to_gene_TSS")]
  
  dtotal$within_assoc_TSS <- FALSE
  dtotal[which(abs(as.numeric(dtotal$diff_to_gene_TSS))<50), "within_assoc_TSS"] <- TRUE
  dtotal$within_assoc_TSS=as.factor(dtotal$within_assoc_TSS)
  
  dme$within_assoc_TSS <- FALSE
  dme[which(abs(as.numeric(dme$diff_to_gene_TSS))<50), "within_assoc_TSS"] <- TRUE
  dme$within_assoc_TSS=as.factor(dme$within_assoc_TSS)
  
  d3$within_assoc_TTS <- FALSE
  d3[which(abs(as.numeric(d3$diff_to_gene_TTS))<50), "within_assoc_TTS"] <- TRUE
  d3$within_assoc_TTS=as.factor(d3$within_assoc_TTS)
  d3=d3[,-which(colnames(d3)=="diff_to_gene_TTS")]
  
  dtotal$within_assoc_TTS <- FALSE
  dtotal[which(abs(as.numeric(dtotal$diff_to_gene_TTS))<50), "within_assoc_TTS"] <- TRUE
  dtotal$within_assoc_TTS=as.factor(dtotal$within_assoc_TTS)
  
  dme$within_assoc_TTS <- FALSE
  dme[which(abs(as.numeric(dme$diff_to_gene_TTS))<50), "within_assoc_TTS"] <- TRUE
  dme$within_assoc_TTS=as.factor(dme$within_assoc_TTS)
  
  # if NA in dist_to_cage_peak, replace it by -11000
  d3[is.na(d3$dist_to_cage_peak),"dist_to_cage_peak"] <- -11000
  dtotal[is.na(dtotal$dist_to_cage_peak),"dist_to_cage_peak"] <- -11000
  
  # if NA in pos_cage_peak, replace it with 0 and if value, replace it with one # added Leo
  d3$pos_cage_peak<-ifelse(is.na(d3$pos_cage_peak),0,1)
  dtotal$pos_cage_peak<-ifelse(is.na(dtotal$pos_cage_peak),0,1)
  
  # if NA in dist_to_polya_site, replace it by -11000
  d3[is.na(d3$dist_to_polya_site),"dist_to_polya_site"] <- -11000
  dtotal[is.na(dtotal$dist_to_polya_site),"dist_to_polya_site"] <- -11000
  
  
  # if NA in predicted_NMD, replace this value by "non_coding"
  d3[is.na(d3$predicted_NMD),"predicted_NMD"] <- "non_coding"
  dtotal[is.na(dtotal$predicted_NMD),"predicted_NMD"] <- "non_coding"
  d3$predicted_NMD=as.factor(d3$predicted_NMD)
  dtotal$predicted_NMD=as.factor(dtotal$predicted_NMD)
  dme[is.na(dme$predicted_NMD),"predicted_NMD"] <- "non_coding"
  dme$predicted_NMD=as.factor(dme$predicted_NMD)
  
  # if NA in Min_sample_cov, replace this value by 0
    d3[is.na(d3$min_sample_cov),"min_sample_cov"] <- 0
  dtotal[is.na(dtotal$min_sample_cov),"min_sample_cov"] <- 0
  
  # if NA in MinCov, replace this value by 0
  d3[is.na(d3$min_cov),"min_cov"] <- 0
  dtotal[is.na(dtotal$min_cov),"min_cov"] <- 0
    
  # if ratioExp is NA, replace this value by 0
  d3[is.na(d3$ratio_exp),"ratio_exp"] <- 0
  dtotal[is.na(dtotal$ratio_exp),"ratio_exp"] <- 0
  
  # if all sdCov is NA, replace by 0. If NA is some of them, replace with median value of sdCov in the dataset
  if (all(is.na(d3$sd_cov))){
    d3$sd_cov = 2
    dtotal$sd_cov = 2
    
  }else{
    medt2=median(as.numeric(d3[!is.na(d3$sd_cov), "sd_cov"]))
    d3[is.na(d3$sd_cov),"sd_cov"] <- medt2
    dtotal[is.na(dtotal$sd_cov),"sd_cov"] <- medt2
    
  }
  
  # if NA in indels, replace this value by 0
  d3[is.na(d3$n_indels),"n_indels"] <- 0
  dtotal[is.na(dtotal$n_indels),"n_indels"] <- 0
  
  # if NA in indels near junctions, replace this value by 0
  d3[is.na(d3$n_indels_junc),"n_indels_junc"] <- 0
  dtotal[is.na(dtotal$n_indels_junc),"n_indels_junc"] <- 0
  
  # Bite
  dtotal[is.na(dtotal$bite), "bite"] <- FALSE
  
  if(opt$exon_cov==0){
    if(!is.null(d1$FIRST_EXON_COV)){
      d3=d3[,-grep('FIRST_EXON_COV',colnames(d3))]
      d3=d3[,-grep('FIRST_EXON_NUMBER',colnames(d3))]
      d3=d3[,-grep('LAST_EXON_NUMBER',colnames(d3))]
      d3=d3[,-grep('LAST_EXON_COV',colnames(d3))]
      d3=d3[,-grep('EXON_min_cov',colnames(d3))]
    }
  }
  if(opt$exon_cov==1){
    if(is.null(d3$FIRST_EXON_COV)|is.null(d3$FIRST_EXON_NUMBER)|is.null(d3$LAST_EXON_COV)|is.null(d3$LAST_EXON_NUMBER)|is.null(d3$EXON_min_cov)){
      print("The information on the exon coverage is not available (columns missing). The variables 'FIRST_EXON_COV','FIRST_EXON_NUMBER','LAST_EXON_NUMBER,'LAST_EXON_COV','EXON_min_cov' will not be considered")
    }
  }
  
  ########################################
  ############# Preprocessing ############
  ########################################
  
  print("-------------------------------------------------")
  print("Preprocessing: remove of the near zero variables")
  
  nzv=nearZeroVar(d3)
  nodel = which(colnames(d3) %in% c("structural_category", "subcategory", "all_canonical"))
  if(length(colnames(d3)[nzv])!=0){
    nzv = setdiff(nzv, nodel)  
  }
    
  if(length(colnames(d3)[nzv])!=0){
    d27=d3[,-nzv]
    print("removed columns: ")
    print(colnames(d3)[nzv])}else{d27=d3
      print("no near zero variances! ")}
  
  print("Preprocessing: Identification of near zero variance predictors")
  dmatrix=d27
  if(!is.null(dmatrix$FSM_class)){
    dmatrix=dmatrix[,-grep('FSM_class',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$structural_category)){
    dmatrix=dmatrix[,-grep('structural_category',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$subcategory)){
    dmatrix=dmatrix[,-grep('subcategory',colnames(dmatrix))]
  }
  if(!is.null(dmatrix$all_canonical)){
    dmatrix=dmatrix[,-grep('all_canonical',colnames(dmatrix))]
  }
  
  
  d=data.matrix(dmatrix)
  r=as.vector(which(apply(d,2,function(x) (anyNA(x)))))
  if (length(r)>0){  d=d[,-r] }
  
  descrCorr=cor(d)
  highCorr <- findCorrelation(descrCorr, cutoff=0.9)
  findCorrelation(descrCorr, cutoff=0.9,verbose = TRUE, names=TRUE)
  
  print("Preprocessing: removing the features with a correlation > 0.9")
  
  d28=d27
  if(length(highCorr)>0){
    print("List of the features removed: ")
    # THE CORRECT PROCEDURE WOULD BE:
    VecBig<-c(colnames(d))
    vec<- VecBig[highCorr]
    for (a in 1: length(vec)){
      print(vec[a])
    }
    d28=d27[,!(colnames(d27) %in% vec)]
  }else{
    print("No feature removed")
  }
  
  
  ###############################
  # training set:               #
  # FSM as positives            #
  # NNC-NC as negatives         #
  ###############################
  
  print("-------------------------------------------------")
  print("Creating positive and negative sets for training and testing")
  
  if(is.null(opt$TP)){
    print("TP: full splice match with support at TSS and TTS (<50 bp distance to annotated one or CAGE/polyA motif support)")
    print("TN: novel not in catalog non canonical")
    fsm=d28[rownames(d3[which(d3$structural_category=="full-splice_match" & 
                                (d3$within_assoc_TSS == TRUE | d3$within_cage_peak == "True" ) & 
                                (d3$within_assoc_TTS == TRUE | d3$found_polyA == "Motif_found" ) ),]),]
                                
    nncnc=d28[rownames(d3[d3$structural_category=="novel_not_in_catalog" & d3$all_canonical=='non_canonical',]),]
    trainingsetcomplet=rbind(fsm,nncnc)
    Class=factor(c(rep("POS",length(fsm$length)),rep("NEG",length(nncnc$length))))
  }else{
    TP=read.table(opt$TP,as.is=TRUE)
    TN=read.table(opt$TN,as.is=TRUE)
    Class=factor(c(rep("POS",length(TP$V1)),rep("NEG",length(TN$V1))))
    trainingset=d28[c(TP$V1,TN$V1),]
  }
  dim(trainingsetcomplet)
  trainingset=trainingsetcomplet[,-grep('structural_category',colnames(trainingsetcomplet))]
  trainingset=trainingset[,-grep('all_canonical',colnames(trainingset))]
  trainingset=trainingset[,-grep('subcategory',colnames(trainingset))]
  trainingset=trainingset[,-grep('ORF_seq',colnames(trainingset))] #Added Leo
  
  #print(summary(trainingset))
  #print(head(trainingset$coding))
  
  
  ###############################
  #########  Partition ##########
  ###############################
  
  print("Partition train set / test set")
  print("Proportion of the label data used for the training (between 0 and 1):")
  print(opt$pourcent_training)
  
  set.seed(3456)
  inTraining=createDataPartition(Class,p=opt$pourcent_training,list=FALSE,times=1)
  
  training=trainingset[inTraining,]
  testing=trainingset[-inTraining,]
  
  print("Description of the training set:")
  print("Table number of positive and negative examples in the training set")
  print(table(trainingsetcomplet[inTraining,]$structural_category))
  print("Table number of positive and negative examples in the testing set")
  print(table(trainingsetcomplet[-inTraining,]$structural_category))
  
  
  
  ####################################
  ##### Machine learning aproach #####
  ####################################
  
  #10 times 10 cross validation
  ctrl=trainControl(method="repeatedcv",repeats=10,
                    classProbs = TRUE,
                    summaryFunction = twoClassSummary,
                    sampling='down',returnData=TRUE,
                    savePredictions = TRUE, returnResamp='all')
  print("-------------------------------------------------")
  print("Creating Random Forest Classifier")
  print("Random Forest parameters:")
  print("-Down sampling in training set")
  print("-10 cross validation")

  set.seed(1)
  #default: 500 trees
  randomforest <- train(training,Class[inTraining],
                        method ='rf',
                        tuneLength=15,
                        metric = "ROC",
                        trControl = ctrl)
  
  imp=varImp(randomforest,scale=FALSE)
  imp<-imp$importance # added Leo
  imp<-data.frame(rownames(imp),imp) # added Leo
  imp<-imp[order(-imp$Overall),] # added Leo
  imp<-imp[,-1, drop=FALSE] # added Leo
  
  write.table(imp,file=paste(opt$dir,'/VarImptable.txt',sep='' ), quote=F, sep="\t") # Modified Leo
  
  
  # Plotting Variable importance
  
  for(i in 1:length(colnames(training))){
    #distribution on the test set
    filename=paste(opt$dir,"/boxplot_",colnames(training[i]),'_training',sep='')
    # if(colnames(training)[i] %in% profile.1$optVariables[1:numbervariable]){
    #   filename=paste(filename,'_selected',sep='')
    # }
    filename=paste(filename,'.pdf',sep='')
    if(is.factor(training[,i])){
      filename=paste(opt$dir,"/hist_",colnames(training[i]),sep='')
      filename=paste(filename,'.pdf',sep='')
      pdf(filename)
      plot(training[,i]~Class[inTraining],main=colnames(training[i]),xlab="test set")
      dev.off()
    }  else{
      if(typeof(training[,i]) == "integer"|typeof(training[,i]) == "double"){
        pdf(filename)
        tmp=training[,i]+1
        boxplot(tmp~Class[inTraining],log='y',main=colnames(training[i]),xlab="test set")
        dev.off()
      }else{
        if(typeof(training[,i])=="logical"){
          filename=paste(opt$dir,"/barplot_",colnames(training[i]),sep='')
          filename=paste(filename,'.pdf',sep='')
          pdf(filename)
          plot(as.factor(training$bite)~Class[inTraining],main=colnames(training[i]),xlab="test set")
          dev.off()
        }else{
          pdf(filename)
          tmp=training[,i]+1
          boxplot(tmp~Class[-inTraining],log='y',main=colnames(training[i]),xlab="test set")
          dev.off()
        }
      }
    }
  }
  
    
  
  ###############################
  ## Classifier in testing set ##
  ###############################
  print("-------------------------------------------------")
  print("Classifier performance in test set")
  test_pred_prob=predict(randomforest,testing,type='prob')
  pred=factor(ifelse(test_pred_prob$POS>=opt$threshold,"POS","NEG"))
  a=data.frame(pred,obs=Class[-inTraining],POS=test_pred_prob$POS,NEG=test_pred_prob$NEG)
  print("AUC, Sens, Spec of the test set")
  print(twoClassSummary(a,lev=levels(a$obs)))
  write("AUC, Sens, Spec of the test set",file=paste(opt$dir,'/statistics_testSet.txt',sep=''))
  write(twoClassSummary(a,lev=levels(a$obs)),file=paste(opt$dir,'/statistics_testSet.txt',sep=''), append=TRUE)
  
  # print("Area under precision-recall curve, Precision, Recall, F ")
  # print(prSummary(a,lev=levels(a$obs)))
  # write("Area under precision-recall curve, Precision, Recall, F ",file=paste(opt$dir,'/statistics_testset.txt',sep=''),append=TRUE)
  # write(prSummary(a,lev=levels(a$obs)),file=paste(opt$dir,'/statistics_testset.txt',sep=''), append=TRUE)
  
  testpredandclass=cbind(test_pred_prob,Class[-inTraining])
  write.table(testpredandclass,paste(opt$dir,'/Pred_test_and_class.txt',sep=''), quote=F, sep="\t")
  
  cm=confusionMatrix(data=pred,reference=Class[-inTraining],positive="POS")
  info="rows:predictions \ncolumns:reference"
  write.table(paste(info,"\n"),file=paste(opt$dir,"/confusionMatrix_testSet.txt",sep=''), quote = F, col.names = F, row.names = F)
  write.table(cm$table,file=paste(opt$dir,"/confusionMatrix_testSet.txt",sep=''), append = TRUE, quote = F, col.names = T, row.names = T)
  write.table("\n",file=paste(opt$dir,"/confusionMatrix_testSet.txt",sep=''),  append = TRUE,quote = F, col.names = F, row.names = F)
  write.table(cm$overall,file=paste(opt$dir,"/confusionMatrix_testSet.txt",sep=''),append=TRUE, quote=F, col.names = F)
  write.table(cm$byClass,file=paste(opt$dir,"/confusionMatrix_testSet.txt",sep=''),append=TRUE, quote=F, col.names = F)
  
  ###############################
  ########## ROC curve ##########
  ###############################
  
  # 1) in function of the probability on the test set (not the same proportion of positives and negatives)
  fileroc=paste(opt$dir,"/ROC_curve_testset.pdf",sep='')
  pdf(fileroc)
  r=roc(as.numeric(Class[-inTraining]),test_pred_prob$POS,percent = TRUE)
  auc(r)
  plot.roc(r)
  text(20,10,paste('AUC =',signif(auc(r),4)))
  text(20,5,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'))
  #dev.off()
  # 2) same proportion positives and negatives on the test set:
  #testing
  #list of the testing positives
  alltestpos=which(Class[-inTraining]=='POS')
  alltestneg=which(Class[-inTraining]=='NEG')
  nbpos=length(alltestpos)
  nbneg=length(alltestneg)
  if(nbpos<nbneg){
    sampleneg=sample(alltestpos,nbpos,replace=FALSE)
    newtest=testing[c(sampleneg,alltestpos),]
    Classnewtest=factor(c(rep('NEG',nbpos),rep('POS',nbpos)))
  }else{
    samplepos=sample(alltestpos,nbneg,replace=FALSE)
    newtest=testing[c(samplepos,alltestneg),]
    Classnewtest=factor(c(rep('POS',nbneg),rep('NEG',nbneg)))
  }
  set.seed(1)
  
  test_pred_prob2=predict(randomforest,newtest,type='prob')
  
  r=roc(as.numeric(Classnewtest),test_pred_prob2$POS,percent = TRUE)
  auc(r)
  #Area under the curve: 
  plot.roc(r)
  lines(r,col='red')
  ci(r)
  text(20,20,paste('AUC (same proportion +/-)=',signif(auc(r),4)),col='red')
  text(20,15,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'),col='red')
  dev.off()
  
  
  ###############################
  ## Classifier in our dataset ##
  ###############################
  print("-------------------------------------------------")
  print("Classifier in our dataset")
  
  preditproba=predict(randomforest,dtotal,type='prob')
  colnames(preditproba) = gsub("NEG","NEG_MLprob", colnames(preditproba))
  colnames(preditproba) = gsub("POS","POS_MLprob", colnames(preditproba))
  print( "Random forest prediction done")
  dtotalandclassprob=cbind(dtotal,preditproba)
  
  mm.neg = dtotalandclassprob[dtotalandclassprob$POS<opt$threshold,]
  mm.neg = mm.neg[-which(mm.neg$structural_category %in% c("full-splice_match","incomplete-splice_match")),"isoform"]
  
  dtotalandclassprob[which(dtotalandclassprob$isoform %in% mm.neg),"ML_classifier"] <- "NEG" 
  dtotalandclassprob[is.na(dtotalandclassprob$ML_classifier),"ML_classifier"] <- "POS"
  
  #print(colnames(dtotalandclassprob))
  #print(colnames(dme))
  
  if (nrow(dme)>0){
    dtotal = rbind(dtotalandclassprob, data.frame(dme, NEG_MLprob=NA, POS_MLprob=NA, ML_classifier="NA"))
    print("ML classifier results (NEG: Isoforms classified as Artifacts)")
    dtotal[dtotal$ML_classifier=="NA","ML_classifier"] <- "1exon"
    print(table(dtotal$ML_classifier))
  }
} else{dtotal=data.frame(d, ML_classifier="NA")}


###############################
### INTRA-PRIMING FILTERING ###
###############################

print("-------------------------------------------------")
print("Intrapriming filtering in our dataset")
#dtotal[,"intra-priming"] = dtotal$perc_A_downstream_TTS > as.numeric(opt$intrapriming) & !(dtotal$structural_category %in% c("full-splice_match","incomplete-splice_match") )
dtotal[,"intra-priming"] = dtotal$perc_A_downstream_TTS > as.numeric(opt$intrapriming) & !(dtotal$structural_category %in% c("full-splice_match") )

print("Intra-priming filtered transcripts:")
print(table(dtotal$`intra-priming`))
  

###############################
######## WRITING RESULTS ######
###############################

dtotal[which(dtotal$ML_classifier=="NEG" | dtotal$`intra-priming`==TRUE),"SQANTI_filter"] <- "Artifact"
dtotal[is.na(dtotal$SQANTI_filter),"SQANTI_filter"] <- "Isoform"
fileres=paste(opt$dir,'/',basename(opt$sqanti_classif),"_filterResults.txt",sep='')
write.table(dtotal,fileres,sep='\t',quote = F, row.names = F)

fileres2=paste(opt$dir,'/',basename(opt$sqanti_classif),"_curatedTranscriptome.txt",sep='')
write.table(dtotal[which(dtotal$SQANTI_filter=="Isoform"),"isoform"],fileres2,sep='\t',quote = F, row.names = F, col.names = F)

print("-------------------------------------------------")
print("SQANTI filter results:")
print(table(dtotal$SQANTI_filter))

print("*** SQANTI filtering finished!")


# END




