% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ADvariantExplorer.R
\name{ADvariantExplorer}
\alias{ADvariantExplorer}
\title{ADvariantExplorer}
\description{
access to GWAS and eQTL data, slighly optimized for Alzheimer's disease
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ADvariantExplorer$new()}}
\item \href{#method-getTargetGene}{\code{ADvariantExplorer$getTargetGene()}}
\item \href{#method-getFilteredGwasTable}{\code{ADvariantExplorer$getFilteredGwasTable()}}
\item \href{#method-getFullGwasTable}{\code{ADvariantExplorer$getFullGwasTable()}}
\item \href{#method-geteQTLSummary}{\code{ADvariantExplorer$geteQTLSummary()}}
\item \href{#method-geteqtlStudyNamesForGroup}{\code{ADvariantExplorer$geteqtlStudyNamesForGroup()}}
\item \href{#method-geteQTLsByLocationAndStudyID}{\code{ADvariantExplorer$geteQTLsByLocationAndStudyID()}}
\item \href{#method-getCellTypes}{\code{ADvariantExplorer$getCellTypes()}}
\item \href{#method-clone}{\code{ADvariantExplorer$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this [R6][R6::R6Class] class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{targetGene}}{character, a candidate gene for this AD GWAS locus}

\item{\code{loc.chrom}}{character, the chromosome of the region to consider}

\item{\code{loc.start}}{numeric, the start of the region to consider}

\item{\code{loc.end}}{numeric, the end of the region to consider}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new instance of ADvariantExplorer
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getTargetGene"></a>}}
\if{latex}{\out{\hypertarget{method-getTargetGene}{}}}
\subsection{Method \code{getTargetGene()}}{
accessor for the object's targetGene field
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$getTargetGene()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the current value of the targetGene member
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getFilteredGwasTable"></a>}}
\if{latex}{\out{\hypertarget{method-getFilteredGwasTable}{}}}
\subsection{Method \code{getFilteredGwasTable()}}{
access to the EMBL-EBI GWAS catalog via bioc, filtered
by location and/or targetGene
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$getFilteredGwasTable(
  targetGeneOnly = TRUE,
  studyString = "alzheimer",
  trim.columns = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{targetGeneOnly}}{logical, default TRUE, otherwise all variants in this regeion,}

\item{\code{studyString}}{character, default "alzheimer", possibly empty string or any other
  study field matchabel substring
affecting any gene whatsoever}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getFullGwasTable"></a>}}
\if{latex}{\out{\hypertarget{method-getFullGwasTable}{}}}
\subsection{Method \code{getFullGwasTable()}}{
access to the complete EMBL-EBI GWAS catalog via bioc
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$getFullGwasTable(trim.columns = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{trim.columns}}{logical reduce to a legible 7 columns}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-geteQTLSummary"></a>}}
\if{latex}{\out{\hypertarget{method-geteQTLSummary}{}}}
\subsection{Method \code{geteQTLSummary()}}{
access to the complete EMBL-EBI eQTL Catalogue
via the catalogueR R package, focusing on expression cis-QTLs
and sQTL (splicing)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$geteQTLSummary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-geteqtlStudyNamesForGroup"></a>}}
\if{latex}{\out{\hypertarget{method-geteqtlStudyNamesForGroup}{}}}
\subsection{Method \code{geteqtlStudyNamesForGroup()}}{
returns the retrieval-ready study names given
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$geteqtlStudyNamesForGroup(groupMatchingString)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{groupMatchingString}}{character, (part of) a qtl_group string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-geteQTLsByLocationAndStudyID"></a>}}
\if{latex}{\out{\hypertarget{method-geteQTLsByLocationAndStudyID}{}}}
\subsection{Method \code{geteQTLsByLocationAndStudyID()}}{
access to the complete EMBL-EBI eQTL Catalogue
via the catalogueR R package, focusing on expression cis-QTLs
and sQTL (splicing
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$geteQTLsByLocationAndStudyID(
  chrom,
  start,
  end,
  studyIDs,
  targetGene.only = TRUE,
  simplify = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{chrom}}{character}

\item{\code{start}}{numeric}

\item{\code{end}}{numeric}

\item{\code{studyIDs}}{character vector, one or more values from eQTL summary unique_id column}

\item{\code{targetGene.only}}{logical, default TRUE, otherwise result has all eQTL genes in the region}

\item{\code{simplify}}{logical, trims columns to just the crucial 4: rsid, pvalue, gene, samples}

\item{\code{method}}{character, either REST or tabix, REST by default}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a data.frame ordered by increasing pvalue.QTL
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getCellTypes"></a>}}
\if{latex}{\out{\hypertarget{method-getCellTypes}{}}}
\subsection{Method \code{getCellTypes()}}{
returns one gene's cell types, or the entire table
from jiang et al 2020, nine cell types across 10 human (and mouse) brain regions:
https://www.sciencedirect.com/science/article/pii/S2589004220309664
scREAD: A Single-Cell RNA-Seq Database for Alzheimer's Disease
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$getCellTypes(geneSymbol = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{geneSymbol}}{character string, either a HUGO gene symbol, or NA}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a character vector or a data.frame (if geneSymbol param is NA)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ADvariantExplorer$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
