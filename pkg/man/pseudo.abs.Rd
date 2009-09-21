\name{pseudo.abs}
\alias{pseudo.abs}

\title{ Selecting pseudo-absences}
\description{
  Collection of methods for building presence-absence (binary) datasets from presence only data. 
  Binary data is needed for running a large number of the species-cimate envelop models available.
  Please read the BIOMOD manual for further explanations.
}
\usage{
pseudo.abs(coor=NULL, status, env = NULL, strategy = "random", distance = 0, nb.points = NULL, add.pres = TRUE, create.dataset=FALSE, species.name = "SpNoName", plot = FALSE, acol = "grey80", pcol = "red")
}

\arguments{
  \item{coor}{ a 2 columns matrix giving the coordinates of the points (the presences recorded and the whole set of potential absences available). It is not needed for the sre strategy. }
  \item{status}{ a vector containing the presence-absence (1-0) information for the coor data. Any point for
   which a "1" is not given will be awarded a zero by default, thus considered an absence }
  \item{strategy}{ the strategy to use for selecting the pseudo-absences (see the examples below): 
   - random: the absences will be taken at random from the whole set of potential absences
   - per: stands for the perimeter to be drawn around the presences as a whole.
   - squares: same as per but the perimeter is drawn individually around each presence. For this
   strategy, information is needed on the distance wanted (distance argument)
   - circles: same as squares but drawing a round perimeter around each presence. For this
   strategy also, information is needed on the distance wanted (distance argument)
   - sre: sites where the environment is considered to be possibly favourable to the species (ac-
   cording to the \code{\link{sre}} model) are unselected as candidate sites for selecting pseudo-absences. For this strategy, the env argument must be given. }
  \item{distance}{ needed for the squares and circles strategies. The unit is the one of the coor data. }
  \item{env}{ a matrix giving information on the environment of the data points as a set of variables (only needed for the sre strategy) }
  \item{nb.points}{ an option for selecting only a limited number of absences at random from the already strategy-selected pseudo-absences (see examples). 
  The default (nb.points=NULL) keeps all the possible absences according to the strategy selected }
  \item{add.pres}{ if True, the output will be an object also containing the presence information (see the value section below)}
  \item{create.dataset} {type True to create a new dataset in your workspace for direct usage. If set to False, only
   an object is returned containing the row indices to take from the original dataset to render a new one where absences
   have been selected according to the settings of the pseudo.abs() function. See sections details and examples for further explanations}
  \item{species.name}{ The dataset will be stored under the name given by this argument, plus the strategy chosen separated by a dot.
   For example, if you give "larix" in this argument and choose the sre strategy, then the output is stored in the new object Dataset.larix.sre } 
  \item{plot}{ an option for plotting the outup dataset of presences and absences obtained }
  \item{acol}{ if plot=T, the colours wanted to plot the absences }
  \item{pcol}{ if plot=T, the colours wanted to plot the presences }
}
\details{
    The pseudo-absences are created by considering any point where the species was not recorded
  and where the environmental conditions are known to cause potential absence. Feeding the models
  with exceeding numbers of absences can significantly disturb the ability of models to discriminate
  meaningful relationships between climate and species distributions. Moreover, running models on
  such heavy databases is incredibly time consuming.
  In addition, some of the chosen absences might unfortunately represent true presences (this is
  particularlyl likely in the case of incomplete samples). The pseudo-absence data therefore gives
  false information for the estimation of the species-climate relationship. Hence, we propose various
  strategies that seek to remove the spurious effects of using uncorrect pseudo-absences before
  running the models.
}
\value{
  A vector is returned containing the rows to take from the original dataset to constitute the new dataset. A dataset
  can also be created if wanted containing the presences and only the absences desired from the original dataset. 
  The advantage of only the vector being returned is that it is a very light object compared to the production of a whole new dataframe. 
}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{sre}} }
\examples{


fulldata <- data.frame(x=rep(1:100, 100), y=rep(1:100, each=100)) 
presences <- data.frame(x=sample(1:100, 10), y=sample(1:100, 10))


#this is our data : recorded presences and a bank of data for which you have no a priori information
plot(fulldata[,1], fulldata[,2])
par(new=TRUE)
plot(presences[,1], presences[,2], col='red', xlim=c(1,100), ylim=c(1,100), pch=19)


pseudo.data <- pseudo.abs(plot=T, coor=rbind(presences, fulldata), status=c(rep(1,10),rep(0,10000)), strategy="circles", distance=8, create.dataset=TRUE) 
#view what is in the dataset produced : absences, and presences (only if add.pres=TRUE)
Dataset.SpNoName.circles.8[1:20,]


#you can plot it yourself (or set plot=TRUE in the function call)
pseudo.data <- pseudo.abs(plot=T, coor=rbind(presences, fulldata), status=c(rep(1,10),rep(0,10000)), strategy="circles", distance=8, create.dataset=TRUE) 
x11()
plot(Dataset.SpNoName.circles.8[,1],Dataset.SpNoName.circles.8[,2]) 


#use the nb.points argument for only selecting a limited number of points
pseudo.data2 <- pseudo.abs(coor=rbind(presences, fulldata), status=c(rep(1,10),rep(0,10000)), strategy="circles", distance=8, create.dataset=TRUE, nb.points=2000, plot=TRUE) 
length(pseudo.data)
length(pseudo.data2)
  
}

\keyword{ models }
\keyword{ dplot }