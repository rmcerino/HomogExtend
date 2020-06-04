\name{control_vh}
\alias{control_vh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control de valores homogeneizados atipicos
}
\description{
Señala aquellos valores que pueden ser atipicos en su entorno. Es necesario revisarlos para decidir si eliminarlos o no.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
control_vh(data, umbral, n_valor_homogeneo, dist_lw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Base de datos con los valores homogeneizados}
  \item{umbral}{Valor entre 0 y 0.5, representa cuanto puede diferir un valor respecto a su entorno}
  \item{n_valor_homogeneo}{String, va entre comillas (""), nombre del valor homogeneizado (si sale de la func_homog es "m2_coef")}
  \item{dist_lw}{Distancia del vecindario}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
