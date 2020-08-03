\name{valor_actualizado}
\alias{valor_actualizado}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Actualización de valores muestrales a un determinado valor de tipo de cambio.
}
\description{
Paso 2 de 4:

Con las elasticidades y la base de datos obtenidas en el paso 1, se procede a la actualización de valores muestrales a un determinado valor de tipo de cambio.
}
\usage{
valor_actualizado(tc_act, datos, elasticidad)
}
\arguments{
  \item{tc_act}{Valor de tipo de cambio al cual se desea actualizar}
  \item{datos}{Base de datos muestral obtenida de la función "elasticidades", es importante que tenga definidos los cuantiles de cada observación}
  \item{elasticidad}{Data frame obtenido de la función "elasticidades", la cual posee las elasticidades por cada cuantil}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
La función arroja la base de datos con dos nuevas variables:

elasticidad = Elasticidad de tipo de cambio en cada observación

valor_actualizado = Es el valor en pesos por metro cuadrado llevado al tipo de cambio de referencia

}

\note{
Este es el segundo paso para actualizar y homogeneizar la base muestral.

Primero: elasticidades()

Segundo: valor_actualizado()

Tercero: parametros_homog()

Cuarto: func_homog()
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
