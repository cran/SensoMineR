\name{format_holos}
\alias{format_holos}

\title{Format Holos data}

\description{Format the data exported from the Holos platform.}

\usage{
format_holos(path.data)
}

\arguments{
\item{path.data}{The path corresponding to the folder containing all Holos data (i.e. one folder per subject with different files: "X_comment.txt", "X_data.txt", "X_last.txt", "X_txt").}
}

\value{A list of 4 objects: \code{IDsubjects}, a dataframe containing the concordance between the names of the subjects as given in the Holos experiment and their ID; \code{datadigit}, a list of S (S = number of subjects) dataframes corresponding to the digit-tracking data; \code{datafinal_coord}, a list of S (S = number of subjects) dataframes corresponding to the final configurations data; and \code{datafinal_verb}, a list of S (S = number of subjects) dataframes corresponding to the final verbalization data.}

\seealso{
\code{\link{analyse_holos}}
}

\examples{
\dontrun{
data <- format_holos(path.data = "C:/MyDirectory/")}
}
