% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvlrt-object-S3-methods.R
\name{plot.pvlrt}
\alias{plot.pvlrt}
\title{Plotting method for a pvlrt object}
\usage{
\method{plot}{pvlrt}(x, type = "bubbleplot", ...)
}
\arguments{
\item{x}{a \code{pvlrt} object; an output of function \code{pvlrt}().}

\item{type}{character string determining the type of plot to show.
Available choices are \code{"bubbleplot"} which calls \link{bubbleplot_pvlrt},
\code{"heatmap"} which calls \link{heatmap_pvlrt}, and
\code{"barplot"} which calls \link{barplot.pvlrt}, with the additional arguments
supplied in ...}

\item{...}{additional arguments passed to heatmap_pvlrt or barplot.pvlrt
depending on \code{type}.}
}
\value{
A \link[ggplot2]{ggplot} object.
}
\description{
Plotting method for a pvlrt object
}
\examples{

# 500 bootstrap iterations (nsim) in the example below
# are for quick demonstration only --
# we recommended setting nsim to 10000 (default) or bigger

test1 <- pvlrt(statin46, nsim = 500)
plot(test1, type = "bubbleplot")
plot(test1, type = "barplot")
plot(test1, type = "heatmap")

}
\seealso{
\link{pvlrt}; \link{bubbleplot_pvlrt}; \link{heatmap_pvlrt}; \link{barplot.pvlrt}
}
