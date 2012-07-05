




#--------------------------------------------------
#--------Loading data and BIOMOD functions--------#
#--------------------------------------------------


#---------Loading from the BIOMOD package---------#

#you might encounter error messages if you miss some of the necessary packages.
#these are ade4, akima, Design, gam, gbm, Hmisc, MASS, mda, nnet, plyr, randomForest, rpart, reshape.
#If it says that the package 'BIOMOD' was not found, then it is not correctly located 
#It should be at R/R-2.9.0(depending on version)/library
library(BIOMOD)


#load data for the examples
data(Sp.Env)  #the environmental variables (Var1 to Var7) and the species distributions (8 species) 
data(CoorXY)  #coordinates for plotting
data(Future1) #future environmental variables for rendering future projections 


#for loading your own data, follow this example.
#set the working directory to where your data is to be able to read it (R->File->Change Directory)
#here we use a text file, but you can also use csv files (see ?read.csv)
MyData <- read.table("my_data.txt", h=T, sep="\t")


#visuallise the dataset
fix(Sp.Env)

#see ?pseudo.abs for generating a dataset with pseudo-absences if only presences are known for a species.


#-----------------------------------
#--------Running the models--------#
#-----------------------------------


#Initialise the datasets. This will create at a dataframe called DataBIOMOD and an object called Biomod.material used
#to store some usefull information when running the code. Make sure you do not delete these files.
#If you have real independant data, a dataset called DataEvalBIOMOD will also be created.
#Change the names of the data by the actual names of your datasets. Change also the indices to correspond to your need.

Initial.State(Response=Sp.Env[,c(11,13)], Explanatory=Sp.Env[,4:10], IndependentResponse=NULL, IndependentExplanatory=NULL)



#This will run the models, the evaluatuation procedures, and determine the treshold maximing the % of presence and absence correctly predicted.
#Select the models you want to run, and the evaluation criteria you want to use.
#Please read the manual for the usage of the different arguments. Pay particular attention to NbRunEval, DataSplit and NbRepPA.


Models(GLM = T, TypeGLM = "quad", Test = "AIC", GBM = T, No.trees = 3000, GAM = T, CTA = T, CV.tree = 50, 
	    ANN = T, CV.ann = 2, SRE = T, Perc025=T, Perc05=F, MDA = T, MARS = T, RF = T, NbRunEval = 2, DataSplit = 80,
      Yweights=NULL, Roc=T, Optimized.Threshold.Roc=T, Kappa=T, TSS=T, KeepPredIndependent = F, VarImport=5,
	    NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)




#--------------------------------------------------
#--------Output analysis and further steps--------#
#--------------------------------------------------

#look at the predictions
load("pred/Pred_Sp281")
dim(Pred_Sp281)
dimnames(Pred_Sp281)[-1]
Pred_Sp281[1:10,,,]



#Plot the response curves for a model. You need to load it first
load("models/Sp281_GAM_full_rep2")
response.plot(Sp281_GAM_full_rep2, Sp.Env[4:10])

load("models/Sp281_ANN_full")
response.plot(Sp281_ANN_full, Sp.Env[4:10])


#Check the evaluation results for each run
Evaluation.results.TSS




#This function produces a dataframe with predictions in binary and/or filtered format for each species. These are stored out of BIOMOD like the initial predictions.
CurrentPred(BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=T, FiltKappa=T, FiltTSS=T) 

#Determine which is the best model for each species according to the evaluation scores.
PredictionBestModel(method='all', Bin.trans = T, Filt.trans = T)

load("pred/BestModelByKappa")
BestModelByKappa


#Some objects are kept in Rs memory but many are stored dircetly on the hard drive.
#Do think about going into the directories created by BIOMOD to see what objects
#have been produced after each step.




#Plot the results. Pred_Sp281 has already been loaded
load("pred/Pred_Sp281_BinKappa")
load("pred/Pred_Sp281_FiltKappa")

par(mfrow=c(1,3))
level.plot(Pred_Sp281[,"GLM",2,1], CoorXY, show.scale=F, title='probabilities')
level.plot(Pred_Sp281_BinKappa[,"GLM",2,1], CoorXY, show.scale=F, title='binary data')
level.plot(Pred_Sp281_FiltKappa[,"GLM",2,1], CoorXY, show.scale=F, title='filtered data')




#--------------------------------------------------------------------------
#--------Making projections in different resolution, space or time--------#
#--------------------------------------------------------------------------


#BIOMOD enables to project the models in other areas, other resolution or other time-slices.
#In this example, one projection will be carried out according to a climate change scenarios (A1.hadCM3).
#The user can select any model he wants. They just however need to have bben run into the Models function.
#Direct transformation of the probability values into presence/absence can be carried out depending on the threshold selected. 


#Like for the initial data, you need to load the future values of the same environmental variables used than for calibrating the models.
MyFutureData <- read.table("FutureData1.txt", h=T, sep="\t")


#we will use the example dataset
Projection(Proj = Future1[,4:10], Proj.name='Future1',
	GLM = T, GBM = T, GAM = T, CTA = T, ANN = T, SRE = T, Perc025=T, Perc05=F, MDA =T, MARS = T, RF = T, 
	BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=F, FiltKappa=F, FiltTSS=F, repetition.models=T)




#-------multiple future scenarios-------#
#if you have an important number of scenarios to run, you may want to do an automated run.
#This is what it could look like.

directory.name <-"Future_files_directory" #this where the scenario are stored as text files.
future.scenarios <- list.files(path=file.name)


for(i in 1:length(future.scenarios)) {

    future.data <- read.table(paste(getwd(), "/", directory.name, "/", future.scenarios[i], sep=""), h=T, sep="\t")

    Projection(Proj = future.data[,3:6], Proj.name = future.scenarios[i],
	     GLM = T, GBM = T, GAM = T, CTA = T, ANN = T, SRE = T, Perc025=T, Perc05=F, MDA =T, MARS = T, RF = T,
	     BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=T, FiltKappa=T, FiltTSS=T, repetition.models=T)

}
#-------


#---------------------------------------------
#--------Analysing projections outputs--------
#---------------------------------------------



load("proj.Future1/Proj_Future1_Sp281")
load("proj.future1/Proj_Future1_Sp277")


x11()
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=F, title='Current Sp281')
level.plot(Proj_Future1_Sp281[,"RF",1,1], CoorXY, show.scale=F, title='future projection Sp281')
level.plot(DataBIOMOD[,9], CoorXY, show.scale=F, title='Current Sp277')
level.plot(Proj_Future1_Sp277[,"RF",1,1], CoorXY, show.scale=F, title='future projection Sp277')




#Isolate the projections of the best models only in a new object
ProjectionBestModel(Proj.name='Future1', Bin.trans=T, Filt.trans=T, method='all')


#This function will produce several consensus methods to render a single projection for each species. See the manual for
#further explanations on the different possibilities and on the information to give in. 
Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=T, binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)




load("proj.Future1/Total_consensus_Future1")
data <- Total_consensus_Future1

par(mfrow=c(2,5))
par(mar=c(0.6,0.6,2,0.6))
level.plot(data[,1,1], CoorXY, show.scale=F)
level.plot(data[,1,2], CoorXY, show.scale=F)
level.plot(data[,1,3], CoorXY, show.scale=F)
level.plot(data[,1,4], CoorXY, show.scale=F)
level.plot(data[,1,5], CoorXY, show.scale=F)
level.plot(data[,2,1], CoorXY, show.scale=F)
level.plot(data[,2,2], CoorXY, show.scale=F)
level.plot(data[,2,3], CoorXY, show.scale=F)
level.plot(data[,2,4], CoorXY, show.scale=F)
level.plot(data[,2,5], CoorXY, show.scale=F)

load("proj.Future1/Total_consensus_Future1_Bin")
data.bin <- Total_consensus_Future1_Bin

x11()
par(mfrow=c(2,5))
par(mar=c(0.6,0.6,2,0.6))
level.plot(data.bin[,1,1], CoorXY, show.scale=F)
level.plot(data.bin[,1,2], CoorXY, show.scale=F)
level.plot(data.bin[,1,3], CoorXY, show.scale=F)
level.plot(data.bin[,1,4], CoorXY, show.scale=F)
level.plot(data.bin[,1,5], CoorXY, show.scale=F)
level.plot(data.bin[,2,1], CoorXY, show.scale=F)
level.plot(data.bin[,2,2], CoorXY, show.scale=F)
level.plot(data.bin[,2,3], CoorXY, show.scale=F)
level.plot(data.bin[,2,4], CoorXY, show.scale=F)
level.plot(data.bin[,2,5], CoorXY, show.scale=F)






#See the manual and the example (?ProbDensFunc) to see how to use this next function
ProbDensFunc()





#---------------------------------
#--------Further analyses--------#
#---------------------------------


#Function to simulate the migration ability of species in a context of climate change. It takes the current prediction, the future prediction, both in presence/absence.
#it also takes the latitutde and longitude of the sites, and the migration distance expressed in degree. Put also the name of the created dataset using quotes.
#Latitude and Longitude should be the coordinates of the sites (vector format)

#we will run a projection with the original dataset to have it in the same
#format as the Future1 projection. We will use the overall mean consensus

Projection(Proj = Sp.Env[,4:10], Proj.name='Current',
  GLM = T, GBM = T, GAM = T, CTA = T, ANN = T, SRE = T, Perc025=T, Perc05=F, MDA =T, MARS = T, 
  RF = T, BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=T, FiltKappa=T, FiltTSS=T, repetition.models=T)

Ensemble.Forecasting(Proj.name= "Current", weight.method='Roc', PCA.median=T, 
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

load("proj.Future1/Total_consensus_Future1")
load("proj.Current/Total_consensus_Current")

Migration(CurrentPred = Total_consensus_Current[,,1], FutureProj = Total_consensus_Future1[,,1],
X=CoorXY[,1], Y=CoorXY[,2], MaxMigr=1, Pred.Save="Future1.Migration")

Future1.Migration[740:760,]







#Function to estimate the absolute/relative number of species loss/gain/turnover by pixel in a climate change context.
#Mind that this needs presence-absence data.

ProjectionBestModel("Current")

load("proj.Future1/Proj_Future1_BestModelByRoc_Bin")
load("proj.Current/Proj_Current_BestModelByRoc_Bin")

Biomod.Turnover(CurrentPred = Proj_Current_BestModelByRoc_Bin[,1,], 
   FutureProj = Proj_Future1_BestModelByRoc_Bin[,1,], Turnover.Save= "Turnover.2050")

Turnover.2050[740:760,]








#Function to estimate the absolute/relative number of pixel loss/gain py species. Also estimate current ad future range size. 
#It takes the current and future distributions of the species (which could have been constrained by Migration) into presence-absence. Put also the name of the created dataset using quotes.


Biomod.RangeSize(CurrentPred = Proj_Current_BestModelByRoc_Bin[,1,], 
   FutureProj = Proj_Future1_BestModelByRoc_Bin[,1,], SpChange.Save="SpChange.2050")

SpChange.2050$Diff.By.Pixel[740:760,]
SpChange.2050$Compt.By.Species














