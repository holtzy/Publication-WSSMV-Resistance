
#Loading packages
SciViewsPackages <- c("svGUI", "svDialogs") 

#install.packages(SciViews)

library(svDialogs)

## Function to compute and compare the AICc of different model for survreg models
#################
AICcsurvival <- function(modelnames) {
  
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,m$loglik[2])
    ParNum<-c(ParNum,sum(m$df))
    numobs<-c(numobs,length(m$linear.predictors))
    effects<-c(effects,unlist(strsplit(as.character(m$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  choiceAIC<-c("View AICc Table", "Save AICc Table","Quit")
  
  exit=0
  while(exit==0){
    resAIC <- dlgList(choiceAIC, title="Please select an option", multiple = FALSE)$res
    
    if (!length(resAIC)||resAIC=="Quit") {
      exit<-1}
    else if (resAIC=="View AICc Table") {
      AICcsorted
      edit(AICcsorted,title="AICc Table")}
    else {
      filename<-paste('Survival_Analyses_', Sys.Date(),sep="")
      AICctablename<-paste("AICcTable_",filename,".csv",sep="")
      #Output to choose the directory and the name of the AICc Table
      AICctablename2<-dlgSave(default=AICctablename,title='Save AICc table to',filters=dlgFilters[c("R","All"),])$res	
      #Save the AICc table
      write.csv(file=AICctablename2, AICcsorted,row.names=F,quote=F,sep=';')
      exit<-1
    }}#End of the while loop
}

#Function to display/save the AICc table
#AICcsurvival(paste("m",0:1,sep=""))


## Function to compute and compare the AICc of different model for survreg models with a Rmd file
#################

AICcsurvivalRmd <- function(modelnames) {
  
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,m$loglik[2])
    ParNum<-c(ParNum,sum(m$df))
    numobs<-c(numobs,length(m$linear.predictors))
    effects<-c(effects,unlist(strsplit(as.character(m$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  ## Display AICc table
  print(AICcsorted)
  filename<-paste('Survival_Analyses_', Sys.Date(),sep="")
  AICctablename<-paste("AICcTable_",filename,".csv",sep="")
  #Save the AICc table
  write.csv(file=AICctablename, AICcsorted,row.names=F,quote=F)
}

#Function to display/save the AICc table
#AICcsurvivalRmd(paste("m",0:1,sep=""))


## Function to compute and compare the AICc of different model for lmer or glmer models
###################

AICcglmer <- function(modelnames) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,as.numeric(logLik(m)))
    ParNum<-c(ParNum,sum(sapply(VarCorr(m),vpars))+length(fixef(m)))
    numobs<-c(numobs,nrow(model.frame(m)))
    effects<-c(effects,unlist(strsplit(as.character(summary(m)$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  choiceAIC<-c("View AICc Table", "Save AICc Table","Quit")
  
  exit=0
  while(exit==0){
    resAIC <- dlgList(choiceAIC, title="Please select an option", multiple = FALSE)$res
    
    if (!length(resAIC)||resAIC=="Quit") {
      exit<-1}
    else if (resAIC=="View AICc Table") {
      AICcsorted
      edit(AICcsorted,title="AICc Table")}
    else {
      filename<-paste('glmer_Analyses_', Sys.Date(),sep="")
      AICctablename<-paste("AICcTable_",filename,".csv",sep="")
      #Output to choose the directory and the name of the AICc Table
      AICctablename2<-dlgSave(default=AICctablename,title='Save AICc table to',filters=dlgFilters[c("R","All"),])$res	
      #Save the AICc table
      write.csv(file=AICctablename2, AICcsorted,row.names=F,quote=F,sep=';')
      exit<-1
    }}#End of the while loop
}



#Function to display/save the AICc table
#AICcglmer(paste("m",0:1,sep=""))

## Function to compute and compare the AICc of different model for lmer or glmer models within a Rmd file
#######################

AICcglmerRmd <- function(modelnames) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,as.numeric(logLik(m)))
    ParNum<-c(ParNum,sum(sapply(VarCorr(m),vpars))+length(fixef(m)))
    numobs<-c(numobs,nrow(model.frame(m)))
    effects<-c(effects,unlist(strsplit(as.character(summary(m)$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  ## Output AICc table
  print(AICcsorted)
  filename<-paste('glmer_Analyses_', Sys.Date(),sep="")
  AICctablename<-paste("AICcTable_",filename,".csv",sep="")
  #Save the AICc table
  write.csv(file=AICctablename, AICcsorted,row.names=F,quote=F)
  
}

#Function to display/save the AICc table
#AICcglmerRmd(paste("m",0:1,sep=""))


library(svDialogs)
## Function to compute and compare the AICc of different model for lm or glm models
AICcglm <- function(modelnames) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,as.numeric(logLik(m)))
    ParNum<-c(ParNum,attr(logLik(m),"df"))
    numobs<-c(numobs,nrow(model.frame(m)))
    effects<-c(effects,unlist(strsplit(as.character(summary(m)$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  choiceAIC<-c("View AICc Table", "Save AICc Table","Quit")
  
  exit=0
  while(exit==0){
    resAIC <- dlgList(choiceAIC, title="Please select an option", multiple = FALSE)$res
    
    if (!length(resAIC)||resAIC=="Quit") {
      exit<-1}
    else if (resAIC=="View AICc Table") {
      AICcsorted
      edit(AICcsorted,title="AICc Table")}
    else {
      filename<-paste('glm_Analyses_', Sys.Date(),sep="")
      AICctablename<-paste("AICcTable_",filename,".csv",sep="")
      #Output to choose the directory and the name of the AICc Table
      AICctablename2<-dlgSave(default=AICctablename,title='Save AICc table to',filters=dlgFilters[c("R","All"),])$res	
      #Save the AICc table
      write.csv(file=AICctablename2, AICcsorted,row.names=F,quote=F,sep=';')
      exit<-1
    }}#End of the while loop
}



#Function to display/save the AICc table
#AICcglm(paste("m",0:1,sep=""))


#lmer
##########################

library(svDialogs)

AICclmer <- function(modelnames) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  modnames<-c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  m<-c()
  effects<-c()
  for (i in 1:length(modelnames)){	
    m<-eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL<-c(LL,as.numeric(logLik(m)))
    ParNum<-c(ParNum,sum(sapply(VarCorr(m),vpars))+length(fixef(m))+1)
    numobs<-c(numobs,nrow(model.frame(m)))
    effects<-c(effects,unlist(strsplit(as.character(summary(m)$call[2]), "~"))[2])
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight,effects)
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  choiceAIC<-c("View AICc Table", "Save AICc Table","Quit")
  
  exit=0
  while(exit==0){
    resAIC <- dlgList(choiceAIC, title="Please select an option", multiple = FALSE)$res
    
    if (!length(resAIC)||resAIC=="Quit") {
      exit<-1}
    else if (resAIC=="View AICc Table") {
      AICcsorted
      edit(AICcsorted,title="AICc Table")}
    else {
      filename<-paste('lmer_Analyses_', Sys.Date(),sep="")
      AICctablename<-paste("AICcTable_",filename,".csv",sep="")
      #Output to choose the directory and the name of the AICc Table			
      AICctablename2<-dlgSave(default=AICctablename,title='Save AICc table to',filters=dlgFilters[c("R","All"),])$res	
      #Save the AICc table
      write.csv(file=AICctablename2, AICcsorted,row.names=F,quote=F,sep=';')
      exit<-1
    }}#End of the while loop
}

#AICclmer(paste("m",0:13,sep=""))

modelnames <- paste("m",0:10,sep="")
i=1
## Function to compute and compare the AICc of different model for ASReml-R models

AICcASReml <- function(modelnames) {
  nummodels <- length(modelnames)
  
  modnames <- c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  fixed <-c()
  Gside <- c()
  Rside <- c()
  for (i in (1:nummodels)){
    
    #print(modelnames[i])
    m <- eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL <- c(LL,m$loglik)
    dfixed <- length(m$vcoeff$fixed)
    dfrand <- length(m$gammas)
    ParNum<-c(ParNum, dfixed+ dfrand)
    numobs<-c(numobs,m$nedf+ dfixed)
    fixed <-c(fixed,as.character(m$fixed.formula)[3])
    ## No correlation between random effects
    if(length(grep("corh",m$call))==0){
      Gside <-c(Gside, as.character(m$random.formula)[2])
      form <- as.character(m$call)
      pattern <- "dat|asreml|G.param|R.param|^F$"
      form <- form[!grepl(pattern,form)]
      test <- as.vector(form[!form%in%c(m$fixed.formula,m$random.formula)])
      Rside <-c(Rside,ifelse(length(test)==0,"NULL",test))
      ## Correlation between random effects
    }else{
      form <- as.character(m$call)
      pattern <- "dat|asreml|G.param|R.param|^F$"
      form <- form[!grepl(pattern,form)]
      form <- as.vector(form[!form%in%c(m$fixed.formula)])
      ## G- and R-side random effects
      if(length(form)==2|length(form)==3){
        Gside <-c(Gside,form[1])
        Rside <-c(Rside,form[2])
        ## Only G-side random effect
      }
      if(length(form)==1){
        Gside <-c(Gside,form)
        Rside <-c(Rside,"NULL")
      }
    }	
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight, fixed, Gside, Rside)
  
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  
  print(AICcsorted)
  
  
  
  choiceAIC<-c("View AICc Table", "Save AICc Table","Quit")
  
  exit=0
  while(exit==0){
    resAIC <- dlgList(choiceAIC, title="Please select an option", multiple = FALSE)$res
    
    if (!length(resAIC)||resAIC=="Quit") {
      exit<-1}
    else if (resAIC=="View AICc Table") {
      AICcsorted
      edit(AICcsorted,title="AICc Table")}
    else {
      filename<-paste('ASReml-R_Analyses_', Sys.Date(),sep="")
      AICctablename<-paste("AICcTable_",filename,".csv",sep="")
      #Output to choose the directory and the name of the AICc Table
      AICctablename2<-dlgSave(default=AICctablename,title='Save AICc table to',filters=dlgFilters[c("R","All"),])$res	
      #Save the AICc table
      write.csv(file=AICctablename2, AICcsorted,row.names=F,quote=F,sep=';')
      exit<-1
    }}#End of the while loop
}



#Function to display/save the AICc table
#AICcASReml(paste("m",0:1,sep=""))



## Function to compute and compare the AICc of different model for ASReml-R models with a Rmd file

AICcASRemlRmd <- function(modelnames) {
  nummodels <- length(modelnames)
  
  modnames <- c()
  ParNum<-c()
  LL<-c()
  numobs<-c()
  fixed <-c()
  Gside <- c()
  Rside <- c()
  for (i in (1:nummodels)){
    #		print(modelnames[i])
    m <- eval(parse(text = modelnames[i]))
    modnames<-c(modnames,modelnames[i])
    LL <- c(LL,m$loglik)
    dfixed <- length(m$vcoeff$fixed)
    dfrand <- length(m$gammas)
    ParNum<-c(ParNum, dfixed+ dfrand)
    numobs<-c(numobs,m$nedf+ dfixed)
    fixed <-c(fixed,as.character(m$fixed.formula)[3])
    ## No correlation between random effects
    if(length(grep("corh",m$call))==0){
      Gside <-c(Gside, as.character(m$random.formula)[2])
      form <- as.character(m$call)
      pattern <- "dat|asreml|G.param|R.param|^F$"
      form <- form[!grepl(pattern,form)]
      test <- as.vector(form[!form%in%c(m$fixed.formula,m$random.formula)])
      Rside <-c(Rside,ifelse(length(test)==0,"NULL",test))
      ## Correlation between random effects
    }else{
      form <- as.character(m$call)
      pattern <- "dat|asreml|G.param|R.param|^F$"
      form <- form[!grepl(pattern,form)]
      form <- as.vector(form[!form%in%c(m$fixed.formula)])
      ## G- and R-side random effects
      if(length(form)==2|length(form)==3){
        Gside <-c(Gside,form[1])
        Rside <-c(Rside,form[2])
        ## Only G-side random effect
      }
      if(length(form)==1){
        Gside <-c(Gside,form)
        Rside <-c(Rside,"NULL")
      }
    }	
  }
  
  LL<-round(LL,3)
  AIC<-round(-2*LL+2*ParNum,3)
  
  AICci<-(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1))
  DeltaAICci<-AICci-min(AICci)
  
  AICc<-round(-2*LL+2*ParNum+(2*ParNum*(ParNum+1))/(numobs-ParNum-1),3)
  DeltaAICc<-round(AICc-min(AICc),3)
  EvRatio<-round(exp(DeltaAICc/2),3)
  weight <- round(exp(-DeltaAICc/2)/sum(exp(-DeltaAICc/2)),3)
  AICtab<-data.frame(modnames,ParNum,LL,AIC,AICc,DeltaAICc,EvRatio,weight, fixed, Gside, Rside)
  
  AICcsorted<-AICtab[order(as.numeric(DeltaAICci)),]
  row.names(AICcsorted)<-NULL
  #AICcsorted
  #edit(AICcsorted,title="AICc Table")
  
  ## Display AICc table
  print(AICcsorted)
  
  AICcsorted
  filename<-paste('ASReml-R_Analyses_', Sys.Date(),sep="")
  AICctablename<-paste("AICcTable_",filename,".csv",sep="")
  #Save the AICc table
  write.csv(file=AICctablename, AICcsorted,row.names=F,quote=F)
}



#Function to display/save the AICc table
#AICcASRemlRmd(paste("m",0:1,sep=""))

