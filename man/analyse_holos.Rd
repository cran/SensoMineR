\name{analyse_holos}
\alias{analyse_holos}
\title{Analyse Holos data}

\description{
Provide numeric tools and graphical tools to analyse Holos data.
}

\usage{
analyse_holos(data, method, axes = c(1, 2), graph = TRUE, export.res = FALSE)
}

\arguments{
\item{data}{A list of object, as returned by the \code{\link{format_holos}} function.}
\item{method}{A string corresponding to the holistic task realized by the subjects during the experiment: "N" for Napping task, "S" for Sorting task, or "SN" for Sorted Napping task.}
\item{axes}{A length 2 vector specifying the components of the factorial analysis to plot. By default, the first two components are plotted.}
\item{graph}{A boolean specifying if the graphical outputs of the factorial analysis should be plotted or not. By default, graph = TRUE.}
\item{export.res}{A boolean specifying if all the graphical outputs should be exported in the working directory or not. By default, export.res = FALSE. NB: If method = "N", setting this argument to TRUE is the only way to access the individual cognitive processes.}
}

\value{
\item{IDsubjects}{A dataframe containing the concordance between the names of the subjects as given in the Holos experiment and their ID.}
\item{summary.task}{The summary of the task realized by the subjects with pieces of information such as the number of steps performed by each subject, the duration of the task, etc.}
\item{res.FA}{The results of the factorial analysis as returned by the \code{\link{MFA}} function of the \code{\link{FactoMineR}} package for Napping data, by the \code{\link{MCA}} function for Sorting data, and by the \code{\link{HMFA}} function for Sorted Napping data.}
\item{datasets}{All the individuals data sets (digit-tracking data of each subject) and panel data sets (merged final configurations and verbalization).}
}

\seealso{
\code{\link{format_holos}}
}

\references{
Le, M.T., Brard, M. & Le, S. (2016). Holos: A collaborative environment for similarity-based holistic approaches. Behavior Research Methods.\cr
Le, M.T., Husson, F. & Le, S. (2014). Digit-tracking: Interpreting the evolution over time of sensory dimensions of an individual product space issued from Napping and sorted Napping. Food Quality and Preference.}

\examples{
\dontrun{
data(videos)

# Example with Napping data
res.N <- analyse_holos(videos, method = "N", export.res = TRUE)
res.N$summary.task$nbstep.time # number of steps and duration of the task for each subject
res.N$summary.task$freq[[1]] # number of times the first subject moved each stimulus during the task
res.N$res.FA # MFA results that can be customized with the plot.MFA function of FactoMineR
res.N$datasets$digitdata[[1]] # digit-tracking data of the first subject
res.N$datasets$finaldata # Napping data (panel level)

# Example with Sorting data
res.S <- analyse_holos(videos, method = "S")
res.S$res.FA # MCA results that can be customized with the plot.MCA function of FactoMineR
res.S$datasets # Sorting data (panel level)
sorting.data <- apply(res.S$datasets, 2, as.factor) ?
res.fast <- fast(sorting.data)
ConsensualWords(res.fast)
}
}

