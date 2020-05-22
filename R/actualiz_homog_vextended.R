elasticidades <- function(data, cuantiles, name_vut, name_tc, name_valor_m2,
                          dist_lw, f_otras_variables, otras_variables){



  # library(DescTools)
  # library(tibble)
  # library(sf)
  # library(spatialreg)
  # library(expss)
  # library(spdep)


  datos <<- data

  names(datos)[names(datos) == name_vut] <<- "vut_vigente"
  names(datos)[names(datos) == name_tc] <<- "tc"
  names(datos)[names(datos) == name_valor_m2] <<- "valor_m2"


  datos$quant <- DescTools::CutQ(datos$vut_vigente, breaks = quantile(datos$vut_vigente,(seq(0, 1, by = 1/cuantiles))))
  datos$quant <<- datos$quant


  cuantiles = 1
  f_otras_variables = T


  if(cuantiles == 1){

    form <- "log(valor_m2) ~ log(tc)"

    } else {
      form <- "log(valor_m2) ~ log(tc) + log(tc):quant"
      }


  if(f_otras_variables == T) {

    form <- paste(form, otras_variables, sep = " + ")

  }

  form <- as.formula(form)
  form <<- form

  print(form)

  ols = lm(form, datos)

  cord <- st_coordinates(datos)
  d <- dnearneigh(cord, 0, dist_lw)
  dlist <- nbdists(d, coordinates(cord))
  idlist <- lapply(dlist, function(x) 1/x)
  lw <- nb2listw(d, glist=idlist ,style="W" , zero.policy = TRUE)

  # Analisis de la dependencia espacial en los residuos y calculo de multiplicadores de Lagrange
  moran <- lm.morantest(lm(form, datos),lw , zero.policy = TRUE) #H0: Independencia espacial
  moran_lm <- lm.LMtests(ols,lw,test="all", zero.policy = T)

  #Eleccion de modelo
  if (moran$p.value > 0.1){

    print("Se estiman las elasticidades mediante modelo lineal")

    b_sig_lm <- data.frame(summary(ols)["coefficients"])
    b_sig_lm <- rownames_to_column(b_sig_lm)
    names(b_sig_lm)[names(b_sig_lm) == "rowname"] <- "vble"
    remove_rownames(b_sig_lm)

    b_total <- cbind(b_sig_lm, b_sig_lm$coefficients.Estimate) #se la vuelvo a pegar para buscar una sola vez, en la col 6

  } else {
    modelo <- ifelse((moran_lm$SARMA$p.value < 0.1), "modelo SAC", ifelse((moran_lm$RLMerr$p.value < 0.1)
                                                                          & (moran_lm$RLMlag$p.value > 0.1), "modelo SEM",
                                                                          "modelo SAR"))

    print(paste("Se estiman las elasticidades mediante", modelo, sep=" "))


    # Regresi?n espacial
    if (modelo == "modelo SAC"){
      regresion <- sacsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }
    if (modelo == "modelo SAR"){
      regresion <- lagsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }
    if (modelo == "modelo SEM"){
      regresion <- errorsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }

    print(summary(regresion))

    ###Significatividad de las coeficientes estimados
    b_sig <- data.frame(summary.sarlm(regresion)["Coef"])
    b_sig <- rownames_to_column(b_sig)
    names(b_sig)[names(b_sig) == "rowname"] <- "vble"
    remove_rownames(b_sig)
    b_sig <- subset(b_sig, vble != "(Intercept)")

    ###Chequeo que sean estadisticamente signicativos y tengan signo positivo
    if(vlookup("log(tc)", b_sig, 5) > 0.1){
      stop(print("ELASTICIDAD NO SIGNIFICATIVA"))
    }else{
      if(vlookup("log(tc)", b_sig, 2) < 0){
        stop(print("ELASTICIDAD NEGATIVA"))}
      else{print("ELASTICIDAD SIGNIFICATIVA Y CON SIGNO ESPERADO - CUANTIL 1")}
    }

    if(cuantiles > 1){
      for (i in 2:cuantiles){
        if(vlookup(paste("log(tc):quantQ", i, sep = ""), b_sig, 5) > 0.1){
          stop(print(paste("ELASTICIDAD NO SIGNIFICATIVA - CUANTIL", i, sep = " ")))
        }else{
          if(vlookup(paste("log(tc):quantQ", i, sep = ""), b_sig, 2) < 0){
            stop(print(paste("ELASTICIDAD NEGATIVA - CUANTIL", i, sep =" ")))}
          else{
            print(paste("ELASTICIDAD SIGNIFICATIVA Y CON SIGNO ESPERADO - CUANTIL", i, sep = " "))
          }
        }
      }
    }

    # Impacto total - solo en modelos SAC y SAR
    if(modelo == "modelo SAC" | modelo == "modelo SAR"){
      impactos <- impacts(regresion, listw=lw)

      a <- data.frame(impactos$total)
      b_total <- cbind(b_sig,a)

    }else{
      b_total <- cbind(b_sig, b_sig$Coef.Estimate) #se la vuelvo a pegar para buscar una sola vez, en la col 6
    }

  }

  ##########################################################################
  ##### Creacion del data frame final
  elasticidad <<- data.frame(q=as.numeric(1:cuantiles))

  if(cuantiles == 1){
    elasticidad$elasticidad[1] <<- vlookup("log(tc)", b_total, 6)
  }else{
    for(i in 2:cuantiles){
      elasticidad$elasticidad[1] <<- vlookup("log(tc)", b_total, 6)
      elasticidad$elasticidad[i] <<- sum(vlookup(paste("log(tc):quantQ", i, sep = ""), b_total, 6),
                                         vlookup("log(tc)", b_total, 6))
    }
  }

  elasticidad$quant <<- paste("Q", elasticidad$q, sep="")
  elasticidad$q <<- NULL

  dir.create("Elasticidad ")

  save(elasticidad, file="Elasticidad/elasticidad.Rda")

  resultado <- list("Elasticidades por cuantil", elasticidad)
  print(resultado)

}

valor_actualizado <- function(tc_act, datos, elasticidad){

  datos <- left_join(datos, elasticidad, by="quant")
  datos <<- datos

  datos$var_tc <- (tc_act/datos$tc) - 1
  datos$var_tc <- datos$var_tc

  datos$valor_actualizado <- (1 + datos$var_tc * datos$elasticidad) * datos$valor_m2
  datos$valor_actualizado <<- datos$valor_actualizado

  print("PROCESO FINALIZADO - RECORDATORIO: GUARDAR LA BASE")

}

func_homog <- function (data, name_sup, name_frente, name_forma, name_ubicacion_cuadra,
                        name_tipodevalor, name_sit_juridica, name_valor_actualizado,

                        f_sup, f_frente, f_forma, f_ubicacion_cuadra, f_tipodevalor,
                        f_sit_juridica, f_otras_variables, otras_variables,

                        dist_lw, p_valor, parcelas) {

  # library(sf)
  # library(tidyverse)
  # library(spdep)
  # library(expss)
  # library(spatialreg)

  datos <<- data

  if (f_sup == T){names(datos)[names(datos) == name_sup] <<- "p_sup"}
  if (f_frente == T) {names(datos)[names(datos) == name_frente] <<- "largo_frente"}
  if (f_forma == T) {names(datos)[names(datos) == name_forma] <<- "forma"}
  if (f_ubicacion_cuadra == T) {names(datos)[names(datos) == name_ubicacion_cuadra] <<- "ubicacion_cuadra"}
  if (f_tipodevalor == T) {names(datos)[names(datos) == name_tipodevalor] <<- "p_tipodevalor"}
  if (f_sit_juridica == T) {names(datos)[names(datos) == name_sit_juridica] <<- "p_sj"}

  names(datos)[names(datos) == name_valor_actualizado] <<- "valor_actualizado"



  if (f_sup == T) {

    form <- paste ("log(valor_actualizado) ~ log(p_sup) +")

  } else {

    form <- paste ("log(valor_actualizado) ~")}

  if (f_frente == T) {

    form <- paste (form, "log(largo_frente) +", sep = " ")}


  if (f_forma == T) {

    form <- paste (form, "forma +", sep = " ")}


  if (f_ubicacion_cuadra == T) {

    form <- paste (form, "ubicacion_cuadra +", sep = " ")}

  if (f_tipodevalor == T) {

    form <- paste (form, "p_tipodevalor +", sep = " ")}

  if (f_sit_juridica == T) {

    form <- paste (form, "p_sj +", sep = " ")}

  if (f_otras_variables == T) {

    form <- paste (form, otras_variables , sep = " ")}


  form <- as.formula(form)
  form <<- form
  print(form)

  ols = lm(form, datos)

  cord <- st_coordinates(datos)
  d <- dnearneigh(cord, 0, dist_lw)
  dlist <- nbdists(d, coordinates(cord))
  idlist <- lapply(dlist, function(x) 1/x)
  lw <<- nb2listw(d, glist=idlist ,style="W" , zero.policy = TRUE)

  moran <- lm.morantest(lm(form, datos),lw , zero.policy = TRUE) #H0: Independencia espacial
  moran_lm <- lm.LMtests(ols,lw,test="all", zero.policy = T)

  if (moran$p.value > 0.1) {

    print("Se estiman las elasticidades mediante modelo lineal")

    b_sig_lm <- data.frame(summary(ols)["coefficients"])
    b_sig_lm <- rownames_to_column(b_sig_lm)
    names(b_sig_lm)[names(b_sig_lm) == "rowname"] <- "vble"
    remove_rownames(b_sig_lm)

    b_total <- cbind(b_sig_lm, b_sig_lm$coefficients.Estimate) #se la vuelvo a pegar para buscar una sola vez, en la col 6
    names(b_total)[names(b_total)=="b_sig_lm$coefficients.Estimate"] <- "a.total"

  } else {

    modelo <- ifelse((moran_lm$SARMA$p.value < 0.1), "modelo SAC",
                     ifelse((moran_lm$RLMerr$p.value < 0.1) & (moran_lm$RLMlag$p.value > 0.1),
                            "modelo SEM", "modelo SAR"))

    print(paste("Se estiman los parámetros para la función de homogeneización mediante un", modelo, sep=" "))

    # Regresión espacial
    if (modelo == "modelo SAC"){
      regresion <- sacsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }
    if (modelo == "modelo SAR"){
      regresion <- lagsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }
    if (modelo == "modelo SEM"){
      regresion <- errorsarlm(ols, data = datos, listw = lw, zero.policy = T, na.action = na.omit)
    }

    if(modelo == "modelo SAC" | modelo == "modelo SAR") {

      a <- impacts(regresion, listw=lw)
      a <- data.frame (a$total)

      b <- data.frame (summary.sarlm(regresion)[["Coef"]])
      b <- rownames_to_column(b)
      names(b)[names(b) == "rowname"] <- "vble"
      remove_rownames(b)
      b <- subset(b, vble != "(Intercept)")
      b_total <- cbind(b,a)

    } else {

      b <- data.frame (summary.sarlm(regresion)[["Coef"]])
      b <- rownames_to_column(b)
      names(b)[names(b) == "rowname"] <- "vble"
      remove_rownames(b)
      b <- subset(b, vble != "(Intercept)")
      b_total <- cbind(b,b$Estimate)
      names(b_total)[names(b_total)=="b$Estimate"] <- "a.total"


    }


  }


  b_sup <- subset(b_total, vble == "log(p_sup)")
  b_sup <- b_sup[,c("vble","a.total", "Pr...z..")]
  names(b_sup)[names(b_sup)=="a.total"] <- "b"
  names(b_sup)[names(b_sup)=="Pr...z.."] <- "p"
  b_sup$b <- ifelse(b_sup$p <= p_valor, b_sup$b, 0)
  b_sup <<- b_sup

  print(paste("El parámetro de la superficie es", round(b_sup$b, digits = 3), sep = " "))

  b_frente <- subset(b_total, vble == "log(largo_frente)")
  b_frente <- b_frente[,c("vble","a.total", "Pr...z..")]
  names(b_frente)[names(b_frente)=="a.total"] <- "b"
  names(b_frente)[names(b_frente)=="Pr...z.."] <- "p"
  b_frente$b <- ifelse(b_frente$p <= p_valor, b_frente$b, 0)
  b_frente <<- b_frente

  print(paste("El parámetro del frente es", round(b_frente$b, digits = 3), sep = " "))

  b_sig <- subset(b_total, vble == "forma1"  | vble == "ubicacion_cuadra1" |vble == "ubicacion_cuadra2" |
                    vble == "ubicacion_cuadra3"| vble == "p_tipodevalor1"| vble == "p_sj1" )

  b_sig <- b_sig[,c("vble","a.total", "Pr...z..")]
  names(b_sig)[names(b_sig)=="a.total"] <- "b"
  names(b_sig)[names(b_sig)=="Pr...z.."] <- "p"
  b_sig$b <- ifelse(b_sig$p < p_valor, b_sig$b, 0)

  b_sig <<- b_sig

  print(paste("El parámetro de la forma es",
              ifelse (is.na(vlookup("forma1", b_sig, 2)) == T, "no existente",
                      round (vlookup("forma1", b_sig, 2), digits = 2)), sep = " "))
  print(paste("El parámetro de la esquina es",
              ifelse (is.na(vlookup("ubicacion_cuadra1", b_sig, 2)) == T, "no existente",
                      round (vlookup("ubicacion_cuadra1", b_sig, 2), digits = 2)), sep = " "))
  print(paste("El parámetro de ubicación en la cuadra interna es",
              ifelse (is.na(vlookup("ubicacion_cuadra2", b_sig, 2)) == T, "no existente",
                      round (vlookup("ubicacion_cuadra2", b_sig, 2), digits = 2)), sep = " "))
  print(paste("El parámetro de ubicación en la cuadra con salida a dos calles es",
              ifelse (is.na(vlookup("ubicacion_cuadra3", b_sig, 2)) == T, "no existente",
                      round (vlookup("ubicacion_cuadra3", b_sig, 2), digits = 2)), sep = " "))
  print(paste("El parámetro del tipo de valor es",
              ifelse (is.na(vlookup("p_tipodevalor1", b_sig, 2)) == T, "no existente",
                      round (vlookup("p_tipodevalor1", b_sig, 2), digits = 2)), sep = " "))
  print(paste("El parámetro de la situación jurídica es",
              ifelse (is.na(vlookup("p_sj1", b_sig, 2)) == T, "no existente",
                      round (vlookup("p_sj1", b_sig, 2), digits = 2)), sep = " "))


  cero_sup <- ifelse (length(b_sup$b) != 0 ,  b_sup$b,  NA)
  cero_frente <- ifelse (length(b_frente$b) != 0, b_frente$b, NA)
  cero_sig <- ifelse (length(b_sig$b) != 0,
                      as.numeric(apply (as.data.frame(b_sig$b), 2, mean)), NA)
  cero <- as.data.frame(rbind(cero_sig, cero_frente, cero_sup))

  if (mean(cero$V1, na.rm = T) == 0 | is.nan(mean(cero$V1, na.rm = T))==T){
    stop(print("LAS VARIABLES NECESARIAS PARA LA FUNCION DE HOMOGENEIZACION SON NO SIGNIFICATIVAS"))}

  else {

    if (is.na(cero_sup) == F) {
      if (cero_sup > 0) stop(print("VARIABLE SUPERFICIE NO TIENE EL SIGNO ESPERADO"))}

    else {

      if (is.na(cero_frente) == F) {
        if (cero_frente < 0) stop( print("VARIABLE LARGO DE FRENTE NO TIENE EL SIGNO ESPERADO"))}

      else{

        if (is.na (cero_sig) == F) {

          if (is.na(vlookup("forma1", b_sig, 2) == F) &
              vlookup("forma1", b_sig, 2) > 0) {
            stop( print("VARIABLE FORMA NO TIENE EL SIGNO ESPERADO"))}


          if (is.na(vlookup("ubicacion_cuadra1", b_sig, 2) == F) &
              vlookup("ubicacion_cuadra1", b_sig, 2) < 0) {
            stop( print ("VARIABLE ESQUINA NO TIENE EL SIGNO ESPERADO"))}


          if (is.na(vlookup("ubicacion_cuadra2", b_sig, 2) == F) &
              vlookup("ubicacion_cuadra2", b_sig, 2) > 0) {
            stop( print ("VARIABLE INTERNO NO TIENE EL SIGNO ESPERADO"))}


          if (is.na(vlookup("p_tipodevalor1", b_sig, 2) == F) &
              vlookup("p_tipodevalor1", b_sig, 2) < 0) {
            stop( print ("VARIABLE TIPO DE VALOR NO TIENE EL SIGNO ESPERADO"))}


          if (is.na(vlookup("p_sj1", b_sig, 2) == F) &
              vlookup("p_sj1", b_sig, 2) > 0) {
            stop( print ("VARIABLE SITUACION JURIDICA NO TIENE EL SIGNO ESPERADO"))}


        }}}}

  print("Seguimos procesando pero vamos bien")

  dir.create("Coefcientes")

  save(b_sup, file = "Coefcientes/b_sup.Rda")
  save(b_frente, file = "Coefcientes/b_frente.Rda")
  save(b_sig, file = "Coefcientes/b_sig.Rda")

  beta <- as.matrix(b_sig[,c("b")])
  matriz_beta <<- beta


  vbles <- data.frame(id = as.numeric(1:length(datos$forma)))
  if (is.na(vlookup("forma1", b_sig, 2))==F) {
    vbles$forma <- as.numeric(as.character(datos$forma))}

  if (is.na(vlookup("ubicacion_cuadra1", b_sig, 2))==F) {
    vbles$ubicacion_cuadra1 <- ifelse (as.numeric(as.character(datos$ubicacion_cuadra))==1,1,0)}

  if (is.na(vlookup("ubicacion_cuadra2", b_sig, 2))==F) {
    vbles$ubicacion_cuadra2 <- ifelse (as.numeric(as.character(datos$ubicacion_cuadra))==2,1,0)}

  if (is.na(vlookup("ubicacion_cuadra3", b_sig, 2))==F) {
    vbles$ubicacion_cuadra3 <- ifelse (as.numeric(as.character(datos$ubicacion_cuadra))==3,1,0)}

  if (is.na(vlookup("p_tipodevalor1", b_sig, 2))==F) {
    vbles$p_tipodevalor <- as.numeric(as.character(datos$p_tipodevalor))}

  if (is.na(vlookup("p_sj1", b_sig, 2))==F) {
    vbles$p_sj <- as.numeric(as.character(datos$p_sj))}

  largo <- as.numeric(length(vbles$id))
  vbles$id <- NULL

  vbles_muestra <<- vbles

  for(i in 1:largo){
    a <- t(as.matrix(as.numeric(vbles [i,])))
    datos$expon[i] <- a %*% beta}


  mediana_sup <<- round(median(parcelas$p_sup))
  save(mediana_sup, file = "Coefcientes/mediana_sup.Rda")

  mediana_frente <<- round(median(parcelas$largo_frente))
  save(mediana_frente, file = "Coefcientes/mediana_frente.Rda")


  datos$coef <- ((datos$p_sup/mediana_sup) ^ ifelse(length(b_sup$b)!=0, b_sup$b, 0)) *
    ((datos$largo_frente/mediana_frente) ^ ifelse(length(b_frente$b)!=0, b_frente$b, 0)) *
    (exp(datos$expon))

  summary(datos$coef)

  datos$coef <- ifelse (datos$coef < 0.2, 0.2,
                        ifelse (datos$coef > 1.5, 1.5,
                                datos$coef))
  datos$coef <<- datos$coef

  coeficientes_muestra <- summary(datos$coef)

  datos$m2_coef <- datos$valor_actualizado/datos$coef

  datos$m2_coef <<- datos$m2_coef

  vh <- summary(datos$valor_actualizado/datos$coef)

  save(datos, file = "Coefcientes/datos_coef.Rda")

  resultados_muestra <- list("RESUMEN DE COEFICIENTES EN LA MUESTRA",coeficientes_muestra, "RESUMEN DE VALORES HOMOGENEIZADOS",vh)


  vbles2 <- data.frame(id = as.numeric(1:length(parcelas$forma)))
  if (is.na(vlookup("forma1", b_sig, 2))==F) {
    vbles2$forma <- as.numeric(as.character(parcelas$forma))}

  if (is.na(vlookup("ubicacion_cuadra1", b_sig, 2))==F) {
    vbles2$ubicacion_cuadra1 <- ifelse (as.numeric(as.character(parcelas$ubicacion_cuadra))==1,1,0)}

  if (is.na(vlookup("ubicacion_cuadra2", b_sig, 2))==F) {
    vbles2$ubicacion_cuadra2 <- ifelse (as.numeric(as.character(parcelas$ubicacion_cuadra))==2,1,0)}

  if (is.na(vlookup("ubicacion_cuadra3", b_sig, 2))==F) {
    vbles2$ubicacion_cuadra3 <- ifelse (as.numeric(as.character(parcelas$ubicacion_cuadra))==3,1,0)}

  if (is.na(vlookup("p_tipodevalor1", b_sig, 2))==F) {
    vbles2$p_tipodevalor <- 0}

  if (is.na(vlookup("p_sj1", b_sig, 2))==F) {
    vbles2$p_sj <- 0}

  largo2 <- as.numeric(length(vbles2$id))
  vbles2$id <- NULL

  vbles_parcelas <<- vbles2


  for(i in 1:largo2){
    a <- t(as.matrix(as.numeric(vbles2 [i,])))
    parcelas$expon[i] <- a %*% beta}


  parcelas$coef <- ((parcelas$p_sup/mediana_sup) ^ ifelse(length(b_sup$b)!=0, b_sup$b, 0)) *
    ((parcelas$largo_frente/mediana_frente) ^ ifelse(length(b_frente$b)!=0, b_frente$b, 0)) *
    (exp(parcelas$expon))


  parcelas$coef <- ifelse (parcelas$coef < 0.2, 0.2,
                           ifelse (parcelas$coef > 1.5, 1.5,
                                   parcelas$coef))

  parcelas$coef <<- parcelas$coef


  save(parcelas, file = "Coefcientes/parcelas_coef.Rda")

  coeficientes_parcelas <- summary(parcelas$coef)

  reaultados_parcelas <- list("RESUMEN COEFICIENTES PARCELAS", coeficientes_parcelas)

  resultados <- list(resultados_muestra, reaultados_parcelas)

  return(resultados)

}





