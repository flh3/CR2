% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sch25.R
\docType{data}
\name{sch25}
\alias{sch25}
\title{Data from 25 schools (based on the NELS dataset)}
\format{
A data frame with 546 rows and 8 variables:
\describe{
\item{schid}{The school identifier (the grouping variable)}
\item{ses}{Student-level socioeconomic status}
\item{byhomewk}{Total amount of time the student spent on homework per week. 1 = None, 2 = Less than one hour, 3 = 1 hour, 4 = 2 hours, 5 = 3 hours, 6 = 4-6 hours, 7 = 7 - 9 hours, 8 = 10 or more}
\item{math}{Mathematics score.}
\item{male}{Dummy coded gender, 1 = male, 0 = female}
\item{minority}{Dummy coded minority status, 1 = yes, 0 = no}
\item{mses}{Aggregated socioeconomic status at the school level}
\item{mhmwk}{Aggregated time spent on homework at the school level}
}
}
\source{
\url{https://nces.ed.gov/pubs92/92030.pdf}
}
\usage{
data(sch25)
}
\description{
For examining the association between amount homework done per week and math outcome.
}
\keyword{datasets}
