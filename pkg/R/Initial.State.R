`Initial.State` <-
function(Response=NULL, Explanatory=NULL, IndependentResponse=NULL, IndependentExplanatory=NULL, sp.name=NULL)
{

  if(exists("DataEvalBIOMOD")) rm(DataEvalBIOMOD, inherits=TRUE)
  if(is.null(sp.name) && is.vector(Response)) stop("you need to give a 'sp.name' if only one species is selected in 'Response'")

  if((!is.null(IndependentResponse) && is.null(IndependentExplanatory)) ||
    (is.null(IndependentResponse) && !is.null(IndependentExplanatory)))
    stop("Independent data should be entered for both response and explanatory variables")

    if(is.null(IndependentResponse)) {
        if(is.vector(Response) || is.double(Response)) {
            Response <- as.data.frame(Response)
            dimnames(Response) <- list(seq(dim(Response)[1]), sp.name)
        }
    }    
    if(!is.null(IndependentResponse)) {
        if(ncol(as.data.frame(Response)) != ncol(as.data.frame(IndependentResponse))) stop("Independent data should have the same number of columns as the data used for calibration")
        
        nb <- 0
        for(i in 1:ncol(IndependentExplanatory)) if(sum(colnames(IndependentExplanatory)[i]==colnames(Explanatory)) == 1) nb <- nb+1
        if(nb != ncol(Explanatory)) stop("The variable names given in IndependentExplanatory do not correspond to the one used for calibrating the models \n") 
        
        #reorder the variables correctly 
        IndependentExplanatory <- IndependentExplanatory[match(colnames(Explanatory),colnames(IndependentExplanatory))]
            
        if(is.vector(Response) || is.double(Response)) {
            Response <- as.data.frame(Response)
            dimnames(Response) <- list(seq(dim(Response)[1]), sp.name)
            IndependentResponse <- as.data.frame(IndependentResponse)
            dimnames(IndependentResponse) <- list(seq(dim(IndependentResponse)[1]), sp.name)
        }
        DataEvalBIOMOD <- cbind(IndependentExplanatory, IndependentResponse)
        assign("DataEvalBIOMOD", DataEvalBIOMOD, pos=1)
        if(dim(na.omit(DataEvalBIOMOD))[1] != dim(DataEvalBIOMOD)[1]) stop("Evaluation data contain NA, some models may not work")
    }
    
    #convert "characters" into "factor"
    VarTypes <- c()
    for (i in 1:dim(Explanatory)[2]){
      if (class(Explanatory[,i]) == "character"){
        Explanatory[,i] = as.factor(Explanatory[,i])
        warning(paste(colnames(Explanatory)[i]," was converted to factor"))
        }
      VarTypes <- c(VarTypes,class(Explanatory[,i]))
      }

    assign("DataBIOMOD", cbind(Explanatory, Response), pos=1)
    Biomod.material <- list()
    Biomod.material[["NbVar"]] <- dim(Explanatory)[2]
    Biomod.material[["VarNames"]] <- colnames(Explanatory)
    Biomod.material[["VarTypes"]] <- VarTypes
    Biomod.material[["NbSpecies"]] <- dim(Response)[2]
    Biomod.material[["species.names"]] <- colnames(Response) 
    Biomod.material[["Independent.data.set"]] <- ifelse(exists("DataEvalBIOMOD"),TRUE,FALSE)
    assign("Biomod.material", Biomod.material, pos=1)
    
    if(dim(na.omit(DataBIOMOD))[1] != dim(DataBIOMOD)[1]) stop("Data contain NA, some models may not work")
    if(Biomod.material[["NbVar"]] <= 2) warning("Only two explanatory variables are selected") 
}

