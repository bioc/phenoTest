\name{getVars2test}
\alias{getVars2test}
\docType{data}
\title{
Get phenotypic variables that were tested.
}
\description{
Returns an object containing the names of the variables that were tested
when the epheno object was created.
Will return an object of class list. Variables of the same type
(categorical, survival, etc) will be in the same slot of the list. The
slot names are the types of the variables.
}
\keyword{datasets}
\examples{
data(epheno)
getVars2test(epheno)
}
\author{
Evarist Planet
}