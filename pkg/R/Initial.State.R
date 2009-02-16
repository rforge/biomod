`Initial.State` <-
function(Response=NULL, Explanatory=NULL, IndependentResponse=NULL, IndependentExplanatory=NULL, sp.name=NULL)
{

  if(exists("DataEvalBIOMOD")) rm(DataEvalBIOMOD, inherits=T)
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
        if(is.vector(Response) || is.double(Response)) {
            Response <- as.data.frame(Response)
            dimnames(Response) <- list(seq(dim(Response)[1]), sp.name)
            IndependentResponse <- as.data.frame(IndependentResponse)
            dimnames(IndependentResponse) <- list(seq(dim(IndependentResponse)[1]), sp.name)
        }
        DataEvalBIOMOD <- cbind(IndependentExplanatory, IndependentResponse)
        assign("DataEvalBIOMOD", DataEvalBIOMOD, pos=1)
    }

    assign("DataBIOMOD", cbind(Explanatory, Response), pos=1)
    Biomod.material <- list()
    Biomod.material[["NbVar"]] <- dim(Explanatory)[2]
    Biomod.material[["VarNames"]] <- colnames(Explanatory)
    Biomod.material[["NbSpecies"]] <- dim(Response)[2]
    Biomod.material[["species.names"]] <- colnames(Response)   
    assign("Biomod.material", Biomod.material, pos=1)
    
    if(Biomod.material[["NbVar"]] <= 2) warning("Only two explanatory variables are selected") 
}

