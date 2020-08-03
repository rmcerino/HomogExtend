\name{parametros_homog}
\alias{parametros_homog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regresion espacial para estimar los parametros de la funcion de homogeneizacion
}
\description{
Paso 3 de 4:

Lo que hace la funcion es estimar los parametros que conformaran la funcion de homogeneizacion. Lo hace a traves de una regresion lineal (por lo general considerando la dependencia espacial).
}
\usage{
parametros_homog(data, name_sup, name_frente, name_forma, name_ubicacion_cuadra,
                  name_tipodevalor, name_sit_juridica, name_valor_actualizado,
                              
                  f_sup, f_frente, f_forma, f_ubicacion_cuadra, f_tipodevalor,
                  f_sit_juridica, f_otras_variables, otras_variables,
                              
                  dist_lw, p_valor, parcelario)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Base de datos muestrales con valores actualizados del paso 2. Es importante que contenga ademas las variables superficie, frente, forma, ubicacion en la cuadra esquina, interno, salida a dos calles, tipo de valor, situacion juridica y el resto de las variables definidas en "otras_variables". Es importante respetar los detalles expresados mas abajo sobre la clase de cada variable.}

  \item{name_sup}{String - debe indicarse entre comillas (""). Nombre de la variable superficie}
  \item{name_frente}{String - debe indicarse entre comillas (""). Nombre de la variable ancho de frente}
  \item{name_forma}{String - debe indicarse entre comillas (""). Nombre de la variable forma}
  \item{name_ubicacion_cuadra}{String - debe indicarse entre comillas (""). Nombre de la variable ubicacion en la cuadra}
  \item{name_tipodevalor}{String - debe indicarse entre comillas (""). Nombre de la variable tipo de valor}
  \item{name_sit_juridica}{String - debe indicarse entre comillas (""). Nombre de la variable situacion juridica}
  \item{name_valor_actualizado}{String - debe indicarse entre comillas (""). Nombre de la variable valor actualizado, que resulta del paso dos}


  \item{f_sup}{Indica si en la formula se incluye la variable superficie. Variable logica, debe indicarse TRUE o FALSE}
  \item{f_frente}{Indica si en la formula se incluye la variable ancho de frente. Variable logica, debe indicarse TRUE o FALSE}
  \item{f_forma}{Indica si en la formula se incluye la variable forma. Variable logica, debe indicarse TRUE o FALSE}
  \item{f_ubicacion_cuadra}{Indica si en la formula se incluye la variable ubicacion en la cuadra. Variable logica, debe indicarse TRUE o FALSE}
  \item{f_tipodevalor}{Indica si en la formula se incluye la variable tipo de valor. Variable logica, debe indicarse TRUE o FALSE}
  \item{f_sit_juridica}{Indica si en la formula se incluye la variable situacion juridica. Variable logica, debe indicarse TRUE o FALSE}

  \item{f_otras_variables}{Indica si en al formula se incluyen otras variables ademas de superficie, ancho de frente, forma, ubicacion en la cuadra, tipo de valor y situacion juridica. Variable logica, debe indicarse TRUE o FALSE}
  \item{otras_variables}{String - debe indicarse entre comillas (""). Debe indicarse como suma el resto de las variables a implementar en la formula que defina los parametros de la funcion de homogeneizacion - excluir las variables superficie, ancho de frente, forma, ubicacion en la cuadra, tipo de valor y situacion juridica. Por ejemplo:
  otras_variables <- "d_alta + perc_edif + ndbi"}

  \item{dist_lw}{Distancia para generar la matriz W de vecinos}
  \item{p_valor}{P_valor referente al nivel de significacion que se utilizara para admitir las variables en las regresiones}

}
\details{
En la base de datos es necesario que se respeten la clase de variables:

Superficie = Es la superficie en metros cuadrados (variable numerica).

Frente = Es el largo (o ancho) del frente en metros (variable numerica).

Forma = Debe ser de clase "factor" donde 0 indique forma regular y 1 irregular.

Ubicacion en la cuadra = Debe ser de clase "factor" donde 0 indique ubicacion en la cuadra medial, 1 esquina, 2 interno, 3 salida a dos o mas calles.

Tipo de Valor = Debe ser de clase "factor" donde 0 indique venta y 1 oferta.

Situacion Juridica = Debe ser de clase "factor" donde 0 indique con escritura y 1 sin escritura.}
\value{
En el directorio especificado se creara una carpeta con la siguiente informacion:

Valores de los parametros superficie, frente, forma, ubicacion en la cuadra esquina, interno, salida a dos calles, tipo de valor y situacion juridica si existiesen. En total son tres archivos: "b_sup", "b_frente", "b_sig".

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Este es el tercer paso para actualizar y homogeneizar la base muestral.

Primero: elasticidades()

Segundo: valor_actualizado()

Tercero: parametros_homog()

Cuarto: func_homog()
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
