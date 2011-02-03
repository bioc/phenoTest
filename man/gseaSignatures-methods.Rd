\name{gseaSignatures-methods}
\docType{methods}
\alias{gseaSignatures-methods}
\alias{gseaSignatures,ANY,character-method}
\alias{gseaSignatures,ANY,GeneSet-method}
\alias{gseaSignatures,epheno,list-method}
\alias{gseaSignatures,matrix,list-method}
\alias{gseaSignatures,numeric,list-method}
\alias{gseaSignatures,matrix,GeneSetCollection-method}
\alias{gseaSignatures,epheno,GeneSetCollection-method}
\alias{gseaSignatures,numeric,GeneSetCollection-method}
\title{Methods for Function gseaSignatures in Package `phenoTest'}
\description{
Methods for function \code{gseaSignatures} in Package `phenoTest'. For
more information read the function's manual.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY", signatures = "character")}}{
Method for signature of class \code{character}.
}

\item{\code{signature(x = "ANY", signatures = "GeneSet")}}{
Method for signature of class \code{character}.
}

\item{\code{signature(x = "epheno", signatures = "list")}}{
Method for an \code{epheno} object and several signatures stored in an
object of class \code{list}.
}

\item{\code{signature(x = "matrix", signatures = "GeneSetCollection")}}{
Method for an \code{matrix} object and several signatures stored in an
object of class \code{GeneSetCollection}.
}

\item{\code{signature(x = "epheno", signatures = "GeneSetCollection")}}{
Method for an \code{epheno} object and several signatures stored in an
object of class \code{GeneSetCollection}.
}

\item{\code{signature(x = "numeric", signatures = "GeneSetCollection")}}{
Method for an \code{numeric} object and several signatures stored in an
object of class \code{GeneSetCollection}.
}

\item{\code{signature(x = "matrix", signatures = "list")}}{
Method for an \code{matrix} object and several signatures stored in an
object of class \code{list}.
}

\item{\code{signature(x = "numeric", signatures = "list")}}{
Method for an \code{numeric} object and several signatures stored in an
object of class \code{list}.
}
}}
\keyword{methods}
