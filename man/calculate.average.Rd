\name{calculate.average}

\alias{calculate.average}

\title{Calcul de moyenne}
\description{Fonction intermediaire}

\usage{
calculate.average(axeACP)
}

\arguments{
  \item{axeACP}{matrix}
}

\details{
 axeACP: 1st col: identificateur, colonne suivante: coord sur axe, produit, juge (avec 0 pour moyenne)
         les premieres lignes correspondent aux moyennes tous juges confondus
         les suivantes aux notes individuelles, classees par juge
cree tab: tableau avec en ligne les juges, 
           en colonne les produits (ou si c'est une AFM les individus moyens ou partiels ou les 2)
           en profondeur les coord sur les differents axes

}

\author{F Husson}

\keyword{internal}
