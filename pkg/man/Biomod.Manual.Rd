\name{Biomod.Manual}
\alias{Biomod.Manual}

\title{ Calling the Manual }
\description{ Directly opens the pdf manual fom the R session
}

\usage{
Biomod.Manual(manual="Biomod Manual"))
}

\arguments{
\item{manual}{ The name (minus the .pdf extension) of the file you want to open }
}

\details{
The updated (and advised) version is \textit{BiomodTutorial.pdf}. Following this file, you will be
able to make your own simulations and use all the BIOMOD tools.\\
The \textit{Biomod Manual.pdf} file containes all informations about models used in BIOMOD.\\
\textbf{NOTE:} it's possible that some examples not run because of old version of BIOMOD functions.\\
}

\author{ Wilfried Thuiller, Bruno Lafourcade, Damien Georges }

\seealso{ \code{\link{help}}}
\examples{

Biomod.Manual()
Biomod.Manual("BiomodTutorial")
}

\keyword{ models }
