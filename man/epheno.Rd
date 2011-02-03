\name{epheno}
\alias{epheno}
\docType{data}
\title{
epheno object.
}
\description{
Object obatined with ExpressionPhenoTest function.
}
\usage{data(epheno)}
\format{
  The format is:
Formal class 'epheno' [package "phenoTest"] with 8 slots
  ..@ p.adjust.method  : chr "none"
  ..@ assayData        :<environment: 0x1050d5a78> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	5 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:5] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	12 obs. of  5 variables:
  .. .. .. ..$ phenoName : Factor w/ 3 levels "lymph.node.status",..: 3 3 3 3 3 1 1 1 2 3 ...
  .. .. .. ..$ phenoClass: Factor w/ 3 levels "categorical",..: 2 2 2 2 2 1 1 1 3 2 ...
  .. .. .. ..$ phenoType : Factor w/ 3 levels "mean","pval",..: 1 1 1 3 3 1 1 3 3 2 ...
  .. .. .. ..$ meanLabel : Factor w/ 5 levels "[45.2,49.2)",..: 1 2 3 NA NA 4 5 NA NA NA ...
  .. .. .. ..$ survTime  : Factor w/ 1 level "Months2Relapse": NA NA NA NA NA NA NA NA 1 NA ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	1000 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            : list()
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  ..@ annotation       : chr "hgu133a"
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	12 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 2 12 0
  .. .. .. ..$ : int [1:3] 2 10 0
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 1 0 0
}
\examples{
data(epheno)
## maybe str(epheno) ; plot(epheno) ...
}
\keyword{datasets}
