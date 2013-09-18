\name{AMOEBA}
\alias{AMOEBA}
\title{A Multidirectional Optimum Ecotope-Based Algorithm}
\description{A function to calculate spatial clusters using the Getis-Ord local statistic (Ord and Getis, 1995). It searches irregular clusters (ecotopes) on a map through boundaries or grid of it.
}
\usage{
AMOEBA(outc, neig, power = 1, cpu = 1)
}
\arguments{
  \item{outc}{
Numeric vector with the study outcome
}
  \item{neig}{
Neighbours object from "spdep" package
}
  \item{power}{
Number of times to apply AMOEBA
}
  \item{cpu}{
Total number of cpu to run faster the algorithm
}
}
\details{
\var{outc} must be the same lenght as the total number of units in \var{neig}.

\var{power} must be a positive integer. This value tell us how many times the AMOEBA is applied. When you apply the algorithm one time, the algorithm is called AMOEBA (Aldstadt and Getis, 2006; Duque et al., 2011), but when you apply more than one time, the algorithm is called AMOEBA+.

\var{cpu} must be a positive integer. Also, be careful and no put a \var{cpu} bigger than you have in your computer because it still will work but it probably will take long time. This value helps the function with his computational time, because it is used in the parallelization of the algorithm (Widener et al, 2012).
}
\value{
Return a vector with the classification of the outcome. The classification is ordered from lowest (low risk) to high (high risk). Also, it has as much \eqn{3 ^ power} levels of risk.
}
\author{
Guillermo Valles.
}
\note{
\var{power} should satisfy that \eqn{3 ^ power < length ( outc )}, because otherwise it is meaningless and the algorithm will fail. For example in the columbus example, the map have length of 49 so you can put a \var{power} as much of 3. More than this value it will fail.

In \var{cpu} you should use more than 1 CPU when your area of study is huge, because otherwise it will increase the computation time in a small area. If you decide to change the value here an advice, put 80\% of your computer's CPU for more efficient performance of the algorithm. For exapmle, if your computer have 8 CPU, you should input "cpu=6".
}

\references{
Ord, J. K. and Getis, A. (1995). Local spatial autocorrelation statistics: Distributional issues and application.

Aldstadt, J. and Getis, A. (2006). Using AMOEBA to create a spatial weights matrix and identify spatial clusters.

Duque, J. C., Aldstadt, J., Velasquez, E., Franco, J. L., and Betancourt, A. (2011). A computationally efficient method for delineating irregularly shaped spatial clusters.

Widener, M. J., Crago, N. C., and Aldstadt, J. (2012). Developing a parallel computational implementation of AMOEBA.
}

\section{Acknowledgments}{This work is included within the project "A longitudinal multilevel analysis of socioeconomic disparities in cardiovascular diseases: 
questioning past evidence with new methodological approaches" supported by a grant from The Swedish Research Council (#D054740, PI: Juan Merlo).}

\examples{
###########################################
########## Example of boundaries ##########
###########################################
require(AMOEBA)

### Columbus OH spatial analysis data set (spdep package)
data(columbus)
require(maptools)
map<-readShapePoly(system.file('etc/shapes/columbus.shp',package='spdep')[1])

### Application of AMOEBA
res<-AMOEBA(columbus$CRIME,col.gal.nb,1,1)

### Plot
color<-rev(rainbow(length(table(res)),start=0,end=2/6))
plot(map,col=color[as.factor(res)])
title('Clustering of crimes at Columbus (OH)')
names<-c('\nLow','\nMedium','\nHigh')
require(plotrix)
color.legend(map@bbox[1,1],map@bbox[2,1]-0.2,map@bbox[1,2],map@bbox[2,1]-0.4,names,color,align='rb')
}

\keyword{ spatial }
\keyword{ cluster }