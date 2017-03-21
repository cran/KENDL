\name{plot.KRKM}
\alias{plot.KRKM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for KRKM objects
}
\description{
Plot the estimated cumulative distribution function obtained by the \code{KRKM} function.
}

\usage{
\method{plot}{KRKM}(x, conf.int = TRUE, lty = 1, col = 1, lwd=1, xlim=NULL, ylim =
NULL, log="x", xlab = NULL, ylab = "CDF", \dots)
}
\arguments{
  \item{x}{
    An object of class \code{KRKM} which is retured by the \code{KRKM} function.
  }
  \item{conf.int}{
    Determines whether confidence intervals will be plotted.  The default value
    is TRUE.
  }

  \item{lty}{
    An integer specifying line type. The default value is 1.
  }
  \item{col}{
    An integer specifying color type. The default value is 1.
  }
  \item{lwd}{
    A numeric value specifying line width. The default value is 1.
  }

   \item{xlim}{
    A vector specifying the lower and upper boundaries for x values.
  }
  \item{ylim}{
    A vector specifying the lower and upper boundaries for y values.
  }
  \item{log}{
  	One of the standard character strings ``x", ``y", or ``xy" to specify logarithmic horizontal and/or vertical axes. The default value is ``x".
  }

  \item{xlab}{
    Label given to the x-axis with NULL as default.
  }
  \item{ylab}{
    Label given to the y-axis with ``CDF" as default.
  }
  \item{\dots}{
    Other arguments allowed for the general \code{plot} function.
  }
}

\author{
Yuchen Yang, Li Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
KRKM
}
\examples{
# simulated data set from a lognormal regression model
 lod <- exp(rnorm(200, -2.79, 1.12))
 e   <- rnorm(200, 0 , 1.21 )
 t   <- exp(-3.05+0.42*lod+e)
 obs <- pmax(t, lod)
 bdl <- t<lod
# KRKM estimator
 fit <- KRKM(obs=obs, bdl=bdl, lod=lod)
 plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

