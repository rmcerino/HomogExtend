\name{elasticidades}
\alias{elasticidades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cálculo de elasticidades del tipo de cambio
}
\description{
Paso 1 de 3:
Estimación de las elasticidades del tipo de cambio, de acuerdo a una determinada cantidad de cuantiles. Que serán luego utilizadas para la actualización de valores muestrales.
}
\usage{
elasticidades(data, parcelario, cuantiles, name_vut, name_tc, name_valor_m2, dist_lw, f_otras_variables, otras_variables)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Base de datos muestrales (tiene que ser tipo sf)}
  \item{parcelario}{Base parcelaria de datos de predicción (debe tener el vut vigente)}
  \item{cuantiles}{Número de cuantiles que se desean estimar (mayor o igual a 1)}
  \item{name_vut}{String - debe indicarse entre comillas (""). Nombre de la variable que reprsenta el último vut vigente existente en la BASE PARCELARIA}
  \item{name_tc}{String - debe indicarse entre comillas (""). Nombre de la variable tipo de cambio}
  \item{name_valor_m2}{String - debe indicarse entre comillas (""). Nombre de la variable muestral valor por metro cuadrado (sin homogeneizar ni actualizar)}
  \item{dist_lw}{Distancia para generar la matriz W de vecinos}
  \item{f_otras_variables}{Indica si en al fórmula se incluyen otras variables además del tipo de cambio. Variable lógica, debe indicarse TRUE o FALSE}
  \item{otras_variables}{String - debe indicarse entre comillas (""). Debe indicarse como suma el resto de las variables a implementar en la fórmula que defina las elasticidades - excluir el tipo de cambio. Por ejemplo:
  otras_variables <- "log(p_sup) + log(frente) + p_tipodevalor"}
}

\details{
Es importante que la variable tipo de cambio sea el tipo de cambio en el momento que se relevó el dato.
}

\value{

Devuelve un data frame llamado elasticidad, donde se presenta la elasticidad estimada para cada cuantil.

Devuelve una base (llamada datos), donde a la base datos original se le ha incorporado las variable quant (cuantil al que pertenece la observación)
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
%%  ~~who you are~~
}

\note{

Este es el primer paso para actualizar y homogeneizar la base muestral:

Primero: elasticidades()

Segundo: valor_actualizado()

Tercero: func_homog()

}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
