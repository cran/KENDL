\name{KRKM}
\alias{KRKM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the exposure distribution in the presence of detection limits
}
\description{
Calculate the exposure distribution based on the kernel reverse Kaplan-Meier(KRKM) estimator proposed by Yang et al. (2017).

}
\usage{
KRKM(obs, bdl, lod,method='formula',b=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
A vector of observations which are the exposure levels if the subject's exposure level can be detected, and the detection limits otherwise.
}
  \item{bdl}{
A logical vector indicating TRUE 1 if the subject's exposure level is below the detection limit, and FALSE otherwise.
}
  \item{lod}{
A vector of detection limits.
}

  \item{method}{
The method of variance estimation. Either `formula' or `bootstrap' can be used. The default value is `formula'.
}

  \item{b}{
The number of bootstrap replicates for the bootstrap-based variance estimation method. The default value is 1000.
}

}
\details{
The KRKM estimator is a nonparametric estimator for the exposure distribution. It is valid regardless of whether the exposure level and detection limit are independent or not. This method requires the data come from simple random sampling and the detection limit is available for each subject. This function provides either formula-based or bootstrap-based variance estimation. The formula-based variance estimation method is computationally faster than the bootstrap-based method, but may underestimate the variance and thus yield poor coverage
probabilities at the points below which there are few observations.
}
\value{
 \item{unique.obs}{
A vector of unique observations in the argument \code{obs}.
}
 \item{cdf}{A vector of estimates for the cumulative distribution function (cdf) at the points given by \code{unique.obs}.}
 \item{se}{A vector of standard error estimates for the cdf estimates}
 \item{lower}{A vector of 95\% lower confidence limits of the cdf values}
 \item{upper}{A vector of 95\% upper confidence limits of the cdf values}
}

\author{
Yuchen Yang, Li Chen
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
plot.KRKM, summary.KRKM
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
}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

