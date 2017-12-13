\name{videos}
\alias{videos}
\docType{data}

\title{
Videos data obtained with Holos
}

\description{Example of Holos data, as obtained with the \code{\link{format_holos}} function.}

\usage{data(videos)}

\format{
A list of 4 objects: \code{name.subjects}, a dataframe containing the concordance between the names of the subjects as given in the Holos experiment and their ID; \code{datadigit}, a list of S (S = number of subjects) dataframes corresponding to the digit-tracking data; \code{datafinal_coord}, a list of S (S = number of subjects) dataframes corresponding to the final configurations data; and \code{datafinal_verb}, a list of S (S = number of subjects) dataframes corresponding to the final verbalization data.
}

\examples{
data(videos)
}

\seealso{
\code{\link{format_holos} \link{analyse_holos}}
}
