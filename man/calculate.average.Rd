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
 axeACP: 1ère col: identificateur, colonne suivante: coord sur axe, produit, juge (avec 0 pour moyenne)
         les première lignes correspondent aux moyennes tous juges confondus
         les suivantes aux notes individuelles, classées par juge
crée tab: tableau avec en ligne les juges, 
           en colonne les produits (ou si c'est une AFM les individus moyens ou partiels ou les 2)
           en profondeur les coord sur les différents axes

}

\author{François Husson}

\keyword{internal}
