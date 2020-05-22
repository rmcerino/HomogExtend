\name{func_homog}
\alias{func_homog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Función de homogeneización}
\description{
Paso 3 de 3:

Cálculo de los coeficientes de homogeneización (a nivel muestral y parcelario) y de los valores homogeneizados por metro cuadrado de la base muestral.
}
\usage{
func_homog(data, name_sup, name_frente, name_forma, name_ubicacion_cuadra, name_tipodevalor, name_sit_juridica, name_valor_actualizado, f_sup, f_frente, f_forma, f_ubicacion_cuadra, f_tipodevalor, f_sit_juridica, f_otras_variables, otras_variables, dist_lw, p_valor, parcelas)
}
\arguments{
  \item{data}{Base de datos muestrales con valores actualizados del paso 2. Es importante que contenga además las variables superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles, tipo de valor, situación jurídica y el resto de las variables definidas en "otras_variables". Es importante respetar los detalles expresados más abajo sobre la clase de cada variable}

  \item{name_sup}{String - debe indicarse entre comillas (""). Nombre de la variable superficie}
  \item{name_frente}{String - debe indicarse entre comillas (""). Nombre de la variable ancho de frente}
  \item{name_forma}{String - debe indicarse entre comillas (""). Nombre de la variable forma}
  \item{name_ubicacion_cuadra}{String - debe indicarse entre comillas (""). Nombre de la variable ubicación en la cuadra}
  \item{name_tipodevalor}{String - debe indicarse entre comillas (""). Nombre de la variable tipo de valor}
  \item{name_sit_juridica}{String - debe indicarse entre comillas (""). Nombre de la variable situación jurídica}
  \item{name_valor_actualizado}{String - debe indicarse entre comillas (""). Nombre de la variable valor actualizado, que resulta del paso dos}


  \item{f_sup}{Indica si en la fórmula se incluye la variable superficie. Variable lógica, debe indicarse TRUE o FALSE}
  \item{f_frente}{Indica si en la fórmula se incluye la variable ancho de frente. Variable lógica, debe indicarse TRUE o FALSE}
  \item{f_forma}{Indica si en la fórmula se incluye la variable forma. Variable lógica, debe indicarse TRUE o FALSE}
  \item{f_ubicacion_cuadra}{Indica si en la fórmula se incluye la variable ubicación en la cuadra. Variable lógica, debe indicarse TRUE o FALSE}
  \item{f_tipodevalor}{Indica si en la fórmula se incluye la variable tipo de valor. Variable lógica, debe indicarse TRUE o FALSE}
  \item{f_sit_juridica}{Indica si en la fórmula se incluye la variable situación jurídica. Variable lógica, debe indicarse TRUE o FALSE}

  \item{f_otras_variables}{Indica si en al fórmula se incluyen otras variables además de superficie, ancho de frente, forma, ubicación en la cuadra, tipo de valor y situación jurídica. Variable lógica, debe indicarse TRUE o FALSE}
  \item{otras_variables}{String - debe indicarse entre comillas (""). Debe indicarse como suma el resto de las variables a implementar en la fórmula que defina los parámetros de la función de homogeneización - excluir las variables superficie, ancho de frente, forma, ubicación en la cuadra, tipo de valor y situación jurídica. Por ejemplo:
  otras_variables <- "d_alta + perc_edif + ndbi"}

  \item{dist_lw}{Distancia para generar la matriz W de vecinos}
  \item{p_valor}{P_valor referente al nivel de significación que se utilizará para admitir las variables en las regresiones}
  \item{parcelas}{Base parcelaria que será utilizada para la predicción. Es importante que contenga las variables superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles}
    }
\details{
En la base de datos es necesario que se respeten la clase de variables:

Superficie = Es la superficie en metros cuadrados (variable numérica)
Ancho de frente = Es el largo (o ancho) del frente en metros (variable numérica)
Forma = Debe ser de clase "factor" donde 0 indique forma regular y 1 irregular
Ubicación en la cuadra = Debe ser de clase "factor" donde 0 indique ubicación en la cuadra medial, 1 esquina, 2 interno, 3 salida a dos o más calles
Tipo de Valor = Debe ser de clase "factor" donde 0 indique venta y 1 oferta
Situación Jurídica = Debe ser de clase "factor" donde 0 indique con escritura y 1 sin escritura

DEBEN EXISTIR TODAS ESTAS VARIABLES EN LA BASE DE DATOS, DE LO CONTRARIO LA FUNCIÓN ARROJARÁ UN ERROR.

}
\value{

En el directorio especificado se creará una carpeta con la siguiente información:

Valores de los parámetros superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles, tipo de valor y situación jurídica si existiesen.

Valor de la mediana de la superficie y el frente del parcelario.

Una base de datos denominada "datos_coef" en .Rda que contendrá los coeficientes de homogeneización y los valores homogeneizados para cada observación muestral.

Una base denominada "parcelas_coef" en .Rda, donde a la base parcelaria se le agregan los coeficientes de homogeneización de cada observación.

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
Tercero: func_homog()
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
