\name{RcppArmadillo-Functions}
\alias{rcpparma_hello_world}
\alias{HuHuCAR_BT_In}
\alias{HuHuCAR_getData}
\title{Set of functions in example RcppArmadillo package}
\description{
  These four functions are created when
  \code{RcppArmadillo.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcpparma_hello_world()
}
\value{
  \code{rcpparma_hello_world()} does not return a value, but displays a
  message to the console.

  \code{rcpparma_outerproduct()} returns a numeric matrix computed as the
  outer (vector) product of \code{x}.

  \code{rcpparma_innerproduct()} returns a double computer as the inner
  (vector) product of \code{x}.

  \code{rcpparma_bothproducts()} returns a list with both the outer and
  inner products.
  
}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Armadillo C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}

