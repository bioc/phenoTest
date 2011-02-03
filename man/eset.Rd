\name{eset}
\alias{eset}
\docType{data}
\title{
Example data.
}
\description{
Example data of class ExpressionSet. 
}
\usage{data(eset)}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ assayData        :<environment: 0x1050d9390> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	7 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:7] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	286 obs. of  7 variables:
  .. .. .. ..$ PID              : int [1:286] 3 5 6 7 8 9 11 14 15 17 ...
  .. .. .. ..$ GEOaccession     : Factor w/ 286 levels "GSM36777","GSM36778",..: 17 20 21 22 24 25 58 59 60 61 ...
  .. .. .. ..$ lymph.node.status: chr [1:286] "negative" "negative" "negative" "negative" ...
  .. .. .. ..$ Months2Relapse   : int [1:286] 101 118 9 106 37 125 109 14 99 137 ...
  .. .. .. ..$ Relapse          : int [1:286] 0 0 1 0 1 0 0 1 0 0 ...
  .. .. .. ..$ ER.Status        : num [1:286] 0 1 0 0 0 1 1 0 1 1 ...
  .. .. .. ..$ BrainRelapse     : int [1:286] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	16 obs. of  3 variables:
  .. .. .. ..$ Column          : chr [1:16] "ID" "GB_ACC" "SPOT_ID" "Species Scientific Name" ...
  .. .. .. ..$ Description     : Factor w/ 15 levels "","A gene symbol, when one is available (from UniGene).",..: 3 5 15 13 12 1 11 1 10 14 ...
  .. .. .. ..$ labelDescription: chr [1:16] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	1000 obs. of  16 variables:
  .. .. .. ..$ ID                              : Factor w/ 22284 levels "1007_s_at","1053_at",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. ..$ GB_ACC                          : Factor w/ 21129 levels "AF052179","AF061832",..: 93 30 95 97 25 24 96 99 28 20 ...
  .. .. .. ..$ SPOT_ID                         : chr [1:1000] NA NA NA NA ...
  .. .. .. ..$ Species.Scientific.Name         : Factor w/ 2 levels "Homo sapiens",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ Annotation.Date                 : Factor w/ 2 levels "Jul 11, 2007",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ Sequence.Type                   : Factor w/ 4 levels "Consensus sequence",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ Sequence.Source                 : Factor w/ 3 levels "Affymetrix Proprietary Database",..: 1 2 1 2 1 2 1 1 2 1 ...
  .. .. .. ..$ Target.Description              : Factor w/ 21363 levels "Consensus includes gb:AI656011 /FEA=EST /DB_XREF=gi:4739990 /DB_XREF=est:tt42e08.x1 /CLONE=IMAGE:2243462 /UG=Hs.116875 KIAA0156"| __truncated__,..: 16 13 18 20 8 7 19 22 11 4 ...
  .. .. .. ..$ Representative.Public.ID        : Factor w/ 21197 levels "AF052179","AF061832",..: 93 30 95 97 25 24 96 99 28 20 ...
  .. .. .. ..$ Gene.Title                      : Factor w/ 14208 levels "ADP-ribosylation factor 1",..: 35 66 46 60 44 97 96 64 26 33 ...
  .. .. .. ..$ Gene.Symbol                     : Factor w/ 13293 levels "ABCF1","ARF1",..: 20 59 40 53 33 96 94 58 15 18 ...
  .. .. .. ..$ ENTREZ_GENE_ID                  : chr [1:1000] "780" "5982" "3310" "7849" ...
  .. .. .. ..$ RefSeq.Transcript.ID            : Factor w/ 13074 levels "NM_000409","NM_000661 /// NM_001024921",..: 37 45 41 52 1 50 49 82 47 4 ...
  .. .. .. ..$ Gene.Ontology.Biological.Process: Factor w/ 7245 levels "","0000074 // regulation of progression through cell cycle // traceable author statement /// 0006139 // nucleobase, nucleoside, nu"| __truncated__,..: 61 22 78 32 79 60 14 63 72 20 ...
  .. .. .. ..$ Gene.Ontology.Cellular.Component: Factor w/ 4148 levels "","0000502 // proteasome complex (sensu Eukaryota) // traceable author statement /// 0005634 // nucleus // inferred from electroni"| __truncated__,..: 72 45 1 44 1 1 42 71 6 68 ...
  .. .. .. ..$ Gene.Ontology.Molecular.Function: Factor w/ 7314 levels "","0000049 // tRNA binding // non-traceable author statement /// 0000166 // nucleotide binding // inferred from electronic annotat"| __truncated__,..: 23 26 27 40 81 18 39 71 74 69 ...
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
  .. .. ..@ data             :'data.frame':	286 obs. of  0 variables
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
\references{
Has been obtained from GEO (GSE2034). Only first 1000 genes where stored
(the rest has been removed).
}
\examples{
data(eset)
## maybe str(eset) ; plot(eset) ...
}
\keyword{datasets}
