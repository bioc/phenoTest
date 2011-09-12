\name{findCopyNumber}
\alias{findCopyNumber}
\title{
  Find copy number regions using expression data in a similar way ACE does.
}
\description{
  Given enrichment scores between two groups of samples and the
  chromosomical positions of those enrichment scores this function finds
  areas where the enrichemnt is bigger/lower than expected if the
  positions where assigned at random.
  Plots of the regions and positions of the enriched regions are provided.
}
\usage{
findCopyNumber(x, minGenes = 15, B = 100, p.adjust.method = "BH",
pvalcutoff = 0.05, mc.cores = 1, useAllPerm = F, genome = "hg19",
chrLengths, sampleGenome = TRUE, useOneChr = FALSE, useIntegrate =
TRUE)
}
\arguments{
  \item{x}{
    An object of class \code{data.frame} with gene or probe identifiers as
    row names and the following columns: es (the enrichment score), chr (the
    chromosome where the gene or probe belong to) and pos (position in the
    chromosome in megabases).
    It can be obtained (from an epheno object) with the function
    getEsPositions.
  }
  \item{minGenes}{
    Minimum number of genes in a row that have to be enriched to mark the
    region as enriched.
    Has to be bigger than 2.
  }
  \item{B}{
    Number of permuations that will be computed to calculate pvalues.
    If \code{useAllPerm} is FALSE this value has to be bigger than 100. 
    If \code{useAllPerm} is TRUE the computations are much more
    expensive, therefore it is not recommended to use a B bigger than 100.
  }
  \item{p.adjust.method}{
    P value adjustment method to be used. p.adjust.methods provides a list
    of available methods.
  }
  \item{pvalcutoff}{
    All genes with an adjusted p value lower than this parameter will be
    considered enriched.
  }
  \item{mc.cores}{
    Number of cores to be used in the computation. If \code{mc.cores} is
    bigger than 1 the \code{multicore} library has to be loaded.
  }
  \item{useAllPerm}{
    If FALSE for each gene only permutations of genes that are in an area
    with similar density (similar number of genes close to them) are used to
    compute pvalues. 
    If TRUE all permutations are used for each gene.

    We recommend to use the option FALSE after having observed that
    the enrichment can depend on the number of genes that are in the area.

    We recommend to use the option TRUE if the positions of the enrichment
    score are equidistant. Take into account that this option is much slower
    and needs less permutations, therefore a smaller \code{B} is preferred.

    See details for more info.
  }
  \item{genome}{
    Genome that will be used to draw cytobands.
  }
  \item{chrLengths}{
    An object of class \code{numeric} containing chromosome names as names.
    This names have to be the same as the ones used in \code{x$chr}
    If missing the last position is used.
  }
  \item{sampleGenome}{
    If positions are sampled over the hole genome (across chromosomes) or
    within each chromosome. This is TRUE by default. 
  }
  \item{useOneChr}{
    Use only one chromosome to build the distribution under the null
    hypothesis that genes/probes are not enriched. By default this is
    FALSE.
    The chromosome that is used is chosen as follows: after removing small
    chromosomes we select the one closest to the median quadratic distance
    to 0.
    Setting this parameter to TRUE decreases processing time.
  }
  \item{useIntegrate}{
    If we want to use \code{integrate} or \code{pnorm} to compute
    pvalues. The first does not assume any distribution for the
    distribution under the null hypothesis, the second assumes it is
    normally distributed.
  }
}
\details{
  Enrichemnt scores can be either log fold changes, log hazard ratios, log
  variabiliy ratios or any other score.
  
  Within each chromosome a smoothed score for each gene is obtained via
  generalized additive models, the smoothing parameter for each chromosome
  being chosen via cross-validation.
  The obtained smoothing parameter of each chromosome will be used in
  permutations.

  We assessed statistical significance by permuting the positions thrue
  the hole genome.
  If \code{useAllPerm} is FALSE for each gene the permutations of genes
  that are in an area with similar density (distance to tenth gene) are
  used to compute pvalues. We observed that genes with similar densities
  tend to have similar smoothed scores.
  If we set 1000 permutations (\code{B}=1000) scores are permuted thrue 
  the hole genome 10 times (1000/100). For each smoothed scored the
  permutations of the 100 smoothed scores with most similar density
  (distance to tenth gene) are used. Therefore each smoothed score will be
  compared to 1000 smoothed scores obtained from permutations.

  If scores are at the same distance in the genome from each other (for
  instance when we have a score every fixed certain bases) the option
  \code{useAllPerm}=TRUE is recommended. In this case every smoothed score
  is compared to all smoothed scores obtained via permutations.
  In this case having 20,000 genes and setting the paramter \code{B=10}
  would mean that the scores are permuted 10 times times thrue the hole
  genome, obtaining 200,000 permuted smoothed scores. Each observed smoothed
  score will be tested against the distribution of the 200,000 permuted
  smoothed scores.

  Only regions with as many genes as told in \code{minGenes} being
  statistically significant (pvalue lower than parameter
  \code{pvalcutoff}) after adjusting pvalues with the method specified in
  \code{p.adjust.method} will be selected as enriched.
}
\value{
  Plots all chromomes and marks the enriched regions.
  Also returns a \code{data.frame} containing the positions of the
  enriched regions. This output can be passed by to the \code{genesInArea}
  function to obtain the names of the genes that are in each region.
}
\author{
  Evarist Planet
}
\seealso{
  getEsPositions, genesInArea
}
\examples{
data(epheno)
phenoNames(epheno)
mypos <- getEsPositions(epheno,'Relapse')
head(mypos)
#regions <- findCopyNumber(mypos)
#head(regions)
}
\keyword{ ~kwd1 }

