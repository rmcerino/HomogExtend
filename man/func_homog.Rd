\name{func_homog}
\alias{func_homog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Funcion de homogeneizacion}
\description{
Paso 4 de 4:

Calculo de los coeficientes de homogeneizacion (a nivel muestral y parcelario) y de los valores homogeneizados por metro cuadrado de la base muestral.
}
\usage{
func_homog(data, b_sup, b_frente, b_sig, cambiar_forma, cambiar_esquina,
           cambiar_interno, parcelario, coef_parcelas,name_sup, name_frente,
           name_forma, name_ubicacion_cuadra)
}
\arguments{
  \item{data}{Base de datos muestrales con valores actualizados del paso 2. Es importante que contenga ademas las variables superficie, frente, forma, ubicacion en la cuadra esquina, interno, salida a dos calles, tipo de valor, situacion juridica y el resto de las variables definidas en "otras_variables". Es importante respetar los detalles expresados mas abajo sobre la clase de cada variable.}

  \item{b_sup}{Data frame obtenido en el paso 3, que contiene el parametro estimado de la superficie y el p valor.}
  \item{b_frente}{Data frame obtenido en el paso 3, que contiene el parametro estimado del frente y el p valor.}
  \item{b_sup}{Data frame obtenido en el paso 3, que contiene los parametros estimado de y el p valor de las variables: forma, ubicacion en la cuadra, tipo de valor y situacion juridica.}

    \item{cambiar_forma}{Variable logica. En caso que se desee cambiar el parametro estimado para la forma indicar TRUE, sino FALSE.}
    \item{cambiar_esquina}{Variable logica. En caso que se desee cambiar el parametro estimado para la esquina indicar TRUE, sino FALSE.}
    \item{cambiar_interno}{Variable logica. En caso que se desee cambiar el parametro estimado para la ubicacion en la cuadra interna indicar TRUE, sino FALSE.}


  \item{parcelas}{Base parcelaria que sera utilizada para la prediccion. Es importante que contenga las variables que se detallan en la formula: superficie, frente, forma, ubicacion en la cuadra esquina, interno, salida a dos calles (tipo de valor y situacion juridica se vuelven cero en toda la base).}
  \item{coef_parcelas}{Indica si se calculan o no los coeficientes en el parcelario. Variable logica, debe indicarse TRUE o FALSE}

  \item{name_sup}{String - debe indicarse entre comillas (""). Nombre de la variable superficie en el parcelario}
  \item{name_frente}{String - debe indicarse entre comillas (""). Nombre de la variable ancho de frente en el parcelario}
  \item{name_forma}{String - debe indicarse entre comillas (""). Nombre de la variable forma en el parcelario}
  \item{name_ubicacion_cuadra}{String - debe indicarse entre comillas (""). Nombre de la variable ubicacion en la cuadra en el parcelario}
      }
\details{
En la base de datos y parcelario es necesario que se respeten la clase de variables:

Superficie = Es la superficie en metros cuadrados (variable numerica).

Frente = Es el largo (o ancho) del frente en metros (variable numerica).

Forma = Debe ser de clase "factor" donde 0 indique forma regular y 1 irregular.

Ubicacion en la cuadra = Debe ser de clase "factor" donde 0 indique ubicacion en la cuadra medial, 1 esquina, 2 interno, 3 salida a dos o mas calles.

Tipo de Valor = Debe ser de clase "factor" donde 0 indique venta y 1 oferta.

Situacion Juridica = Debe ser de clase "factor" donde 0 indique con escritura y 1 sin escritura.
}
\value{
Dentro de la carpete "Coeficientes" creada en el paso anterior pasara lo siguiente:

En el caso que se cambe algun parametro, se pisara la base b_sig para que el data frame se actualice y se incorporen los cambios necesarios.

Valor de la mediana de la superficie y el frente del parcelario.

Una base de datos denominada "datos_coef" en .Rda que contendra los coeficientes de homogeneizacion y los valores homogeneizados para cada observacion muestral.

Una base denominada "parcelas_coef" en .Rda, donde a la base parcelaria se le agregan los coeficientes de homogeneizacion de cada observacion, en caso de ser aplicable.

Un txt denominado "form_coef" donde se encuentra la funcion de homogeneizacion para el parcelario. Este txt debe ser enviado a sistemas para implementar los coeficientes en el parcelario. Este archivi "form_coef" tambien se guarrda en Rda.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Este es el cuarto paso para actualizar y homogeneizar la base muestral.

Primero: elasticidades()

Segundo: valor_actualizado()

Tercero: parametros_homog()

Cuarto: func_homog()
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
