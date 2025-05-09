% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/methods-shennong.R
\name{sn_layer_data}
\alias{sn_layer_data}
\alias{sn_layer_data<-}
\alias{sn_layer_data,Shennong-method}
\alias{sn_layer_data<-,Shennong-method}
\title{Get Layer Data}
\usage{
sn_layer_data(object, layer = NULL, assay = NULL, ...)

sn_layer_data(object, layer, assay = NULL, ...) <- value

\S4method{sn_layer_data}{Shennong}(object, layer = NULL, assay = NULL, ...)

\S4method{sn_layer_data}{Shennong}(object, layer = NULL, assay = NULL, ...) <- value
}
\arguments{
\item{object}{A Shennong object.}

\item{layer}{Layer name or index.}

\item{assay}{Optional: Assay name.}

\item{...}{Additional arguments.}

\item{value}{The data matrix to assign.}
}
\description{
Get Layer Data

Set Layer Data<-
}
