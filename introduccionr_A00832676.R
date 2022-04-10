# Axel Quiroga Caldera
# A00832676

# Juan Angel Lucio Rojas
# A00833112

# Ejercicio 1
# install.packages("stringi")
# install.packages("random")


library(random)
library(stringi)

N <- sample(1:100,1,replace=F) # Generando un numero aleatorio para la secuencia

generador_ADN <- function(tamanio) {
  secuencia <-  stri_rand_strings(tamanio, 1, '[A C G T]')
  secuencia
  return(secuencia)
}
secuencia <- generador_ADN(tamanio = N) #Llamando la funcion PRUEBA EJERCICIO 1

combined <- paste0(secuencia,collapse = '') #Juntando la cadena  

lista <- as.list(strsplit(combined, "")[[1]]) #Haciendolo lista


# Ejercicio 2
tamanio_ADN <- function(lista) {
  tamanio <- length(lista)
  return(tamanio)
  tamanio
}
tamanio <- tamanio_ADN(lista = secuencia) #Llamando la funcion PRUEBA EJERCICIO 2

# Ejercicio 3
vector_pctg <- c()
porcentaje <- function(cadena_ADN) {
  A_cont <- 0
  T_cont <- 0
  G_cont <- 0
  C_cont <- 0
  n <- length(cadena_ADN)
  for (i in 1:n){
    if (secuencia[i] == "A"){
      A_cont <- A_cont + 1
    }
    else if(secuencia[i] == "T"){
      T_cont <- T_cont + 1
    }
    else if(secuencia[i] == "G"){
      G_cont <- G_cont + 1
    }
    else if(secuencia[i] == "C"){
      C_cont <- C_cont + 1
    }
  }
  
  A_pctg <- (A_cont / n) * 100
  T_pctg <- (T_cont / n) * 100
  G_pctg <- (G_cont / n) * 100
  C_pctg <- (C_cont / n) * 100
  
  vector_pctg <- c(A_pctg, T_pctg, G_pctg, C_pctg)
  sprintf("Adenina: %s Citosina: %s Guanina: %s Timina: %s", vector_pctg[1], vector_pctg[2], vector_pctg[3], vector_pctg[4])
}
porcentaje(cadena_ADN = secuencia) # PRUEBA EJERCICIO 3

# Ejercicio 4
m_arn <- c()

ARN <- function(ADN) {
  molecula_split <- as.list(strsplit(ADN, "")[[1]])
  for(i in 1:length(molecula_split)) {
    if(molecula_split[i] == "A") {
      m_arn <- c(m_arn, "U")
    } 
    else if(molecula_split[i] == "C") {
      m_arn <- c(m_arn, "G")
    } 
    else if(molecula_split[i] == "T") {
      m_arn <- c(m_arn, "A")
    } 
    else if(molecula_split[i] == "G") {
      m_arn <- c(m_arn, "C")
      
    }
  }
  arn_combinada <- paste0(m_arn,collapse = '')
  return(arn_combinada)
  
}
secuencia_ARN <- ARN(ADN = combined) # PRUEBA EJERCICIO 4

# Ejercicio 5
# Utilizando bioseq
# Fuente: https://search.r-project.org/CRAN/refmans/bioseq/html/seq_translate.html
library(bioseq)
secuencia_adn <- generador_ADN(tamanio = 50) #Llamando la funcion de los adns
combined_adn <- paste0(secuencia_adn,collapse = '') # juntandola con paste0
ARN <- ARN(ADN=combined_adn) # hacinedola ARN

aminoacidos <- function(cadena_arn) {
  cadena <- rna(cadena_arn)
  cadena_aminoacidos <- seq_translate(cadena)
  return(cadena_aminoacidos)
}
aminoacidos(cadena_arn = ARN) # PRUEBA EJERCICIO 5

# Ejercicio 6
sentido <- c("5'", "3'")
molecula <- "TGCGATAC"
m_inversa  <- c()
s_inversa <- c()

inversa <- function(ADN, direccion) {
  molecula_split <- as.list(strsplit(ADN, "")[[1]])
  for(i in length(molecula_split):1) {
    m_inversa <- c(m_inversa, molecula_split[i])
  }
  inversa_combinada <- paste0(m_inversa,collapse = '')
  s_inverso <- c(direccion[2], direccion[1])
  
  sprintf("Hebra inversa: %s-%s-%s", s_inverso[1], inversa_combinada, s_inverso[2])
  
}
sprintf("Hebra directa: %s-%s-%s", sentido[1], molecula, sentido[2])
inversa(ADN = molecula, direccion = sentido) # PRUEBA EJERCICIO 6

# Ejercicio 7
sentido_2 <- c("5'", "3'")
molecula_2 <- "TGCGATAC"
m_complementaria <- c()
s_complementaria <- c()

complementaria <- function(ADN, direccion) {
  molecula_split <- as.list(strsplit(ADN, "")[[1]])
  for(i in 1:length(molecula_split)) {
    if(molecula_split[i] == "A") {
      m_complementaria <- c(m_complementaria, "T")
    } 
    else if(molecula_split[i] == "C") {
      m_complementaria <- c(m_complementaria, "G")
    } 
    else if(molecula_split[i] == "T") {
      m_complementaria <- c(m_complementaria, "A")
    } 
    else if(molecula_split[i] == "G") {
      m_complementaria <- c(m_complementaria, "C")
      
    }
  }
  complementaria_combinada <- paste0(m_complementaria,collapse = '')
  s_complementaria <- c(direccion[2], direccion[1])
  
  sprintf("Hebra complementaria: %s-%s-%s", s_complementaria[1], complementaria_combinada, s_complementaria[2])
}

sprintf("Hebra directa: %s-%s-%s", sentido_2[1], molecula_2, sentido_2[2])
complementaria(ADN = molecula_2, direccion = sentido_2) # PRUEBA EJERCICIO 7

# Ejercicio 8
sentido_complementario <- c("3'", "5'")
molecula_complementaria <- "ACGCTATG"
m_com_inversa  <- c()
s_com_inversa <- c()

inversa_complementaria <- function(ADN_complementaria, direccion_complementaria) {
  molecula_split <- as.list(strsplit(ADN_complementaria, "")[[1]])
  for(i in length(molecula_split):1) {
    m_com_inversa <- c(m_com_inversa, molecula_split[i])
  }
  inversa_combinada <- paste0(m_com_inversa,collapse = '')
  s_com_inverso <- c(direccion_complementaria[2], direccion_complementaria[1])
  
  sprintf("Hebra complementaria-inversa: %s-%s-%s", s_com_inverso[1], inversa_combinada, s_com_inverso[2])
  
}
sprintf("Hebra complementaria: %s-%s-%s", sentido_complementario[1], molecula_complementaria, sentido_complementario[2])
inversa_complementaria(ADN_complementaria = molecula_complementaria, direccion_complementaria = sentido_complementario) # PRUEBA EJERCICIO 8

# Ejercicio 9 Pruebas, las pruebas estan abajo de cada funcion :)

