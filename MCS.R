
# matrix_file <- list.files(path = "/Users/benny/MCS/MatriciMM/", pattern = "\\.mat$", full.names = TRUE)

# Funzione per caricare la matrice dal file .mat e restituire la matrice sparsa
library(Matrix)
library(spam)
library(spam64)
library(tictoc)
library(matrixcalc)
library(pryr)
library(here)
library(SparseM)



# Funzione per caricare la matrice dal file .mtx e restituire la matrice sparsa
load_matrix_from_file <- function(filepath) {
  A <- readMM(filepath)
  
  # Verifica se la matrice è sparsa
  if (is.sparseMatrix(A)) {
    cat("La matrice A è sparsa \n")
  }
  return(A)
}

is.sparseMatrix <- function(x) is(x, 'sparseMatrix')

# Definisci il percorso del file .mtx
# filepath <- here("MatriciMM", "ex15.mtx")
filepath <- here ("/Users/benny/MCS/MatriciMM/ex15.mtx")

# Carica e stampa la matrice dal file .mtx
# A <- load_matrix_from_file(filepath)

# Controlla se la matrice è simmetrica controllando i valori non nulli della trasposta

is_symmetric <- function(A) {
  A_transpose <- t(A)
  return(all(A@x == A_transpose@x))
}

#########DA CONTROLLARE ---> bisogna passare un arigomento e renderla funzione, tipo con un booleano


# Crea il vettore 'b' a partire dai valori contenuti nella matrice 'A' passata come parametro affinché la soluzione del
# sistema lineare sia un vettore 'x' composto da tutti 1
create_b_vector <- function(A) {
  n <- nrow(A)
  ones_vector <- rep(1, n)
  b <- A %*% ones_vector
  #DEBUG: print(b)
  return(b)
}

cholesky_decomposition <- function(A) {
  # Ottieni la memoria usata prima della fattorizzazione di Cholesky    
  memoria_iniziale <- mem_used() #usare memory.size() con windows
  # DEBUG cat("MEMORIA INIZIALE CHOL", memoria_iniziale)
  # Calcola la fattorizzazione di Cholesky 
  tryCatch({
    factor <- Cholesky(A)
    cat("La matrice A è definita positiva \n")
  }, error = function(err) {
    cat("La matrice A non è definita positiva \n")
  })
  # Ottieni la memoria usata dopo la fattorizzazione di Cholesky
  memoria_finale <- mem_used() #usare memory.size() con windows
  # DEBUG cat("MEMORIA FINALE CHOL", memoria_finale)
  # Calcola la memoria utilizzata dalla funzione
  memoria_utilizzata_chol <- memoria_finale - memoria_iniziale
  # DEBUG cat("MEMORIA TOTALE CHOL", memoria_utilizzata_chol)
  return(list(factor = factor, memoria_utilizzata_chol = memoria_utilizzata_chol))
}

# Esegui la fattorizzazione di Cholesky della matrice 'A' e salva il risultato in 'result'
# result <- cholesky_decomposition(A)

#########DA CONTROLLARE ---> bisogna passare un arigomento e renderla funzione, tipo con un booleano
# Controlla se 'result' contiene la matrice triangolare inferiore 'L'
# if (!is.null(result$factor)) {
# Stampa la matrice triangolare inferiore 'L'
# print(result$factor)
# } else {
# print("La matrice A non è definita positiva")
# }

# print(result$factor) # Matrice triangolare inferiore L (fattorizzazione di Cholesky)
# print(result$memoria_utilizzata_chol) # Memoria utilizzata durante la fattorizzazione

# print(result$factor) # Matrice triangolare inferiore L (fattorizzazione di Cholesky)
# print(result$memoria_utilizzata_chol) # Memoria utilizzata durante la fattorizzazione
# print(A)

# Risolve il sistema lineare Ax=b prendendo in ingresso la matrice 'A' fattorizzata con il metodo di Cholesky e il
# vettore 'b'. Calcola anche la memoria utilizzata durante il processo
solve_linear_system <- function(factor, b) {
  # Ottieni la memoria usata prima della fattorizzazione di Cholesky
  memoria_iniziale <- mem_used()
  # DEBUG: cat("memoria iniziale", memoria_iniziale)
  x <- solve(factor, b)
  
  # DEBUG: print(x)
  
  # Ottieni la memoria usata dopo la fattorizzazione di Cholesky
  memoria_finale <- mem_used()
  # DEBUG cat("memoria finale", memoria_finale)
  cat("\n\nSistema Risolto\n\n")
  
  # Calcola la memoria utilizzata dalla funzione
  memoria_utilizzata_sistemaLin <- memoria_finale - memoria_iniziale
  # DEBUG: cat("Memoria finale del sistema lineare:", memoria_utilizzata_sistemaLin, "\n")
  return(list(x = x, memoria_utilizzata_sistemaLin = memoria_utilizzata_sistemaLin))
}


compute_relative_error <- function(x) {
  n <- length(x)
  x_esatto <- rep(1, n)
  errore_relativo <- norm(x - x_esatto, type = c("2")) / norm(x_esatto, type = c("2"))
  return(errore_relativo)
}

compute_percent_zeros <- function(A) {
  n_nonzero <- sum(A != 0)
  n_total <- prod(dim(A))
  percent_zero <- 100 * (1 - (n_nonzero / n_total))
  return(percent_zero)
}


# Calcola il numero di elementi diversi da 0 presenti nella matrice A
compute_num_nonzeros <- function(A) {
  nonzeros <- sum(A != 0)
  return(nonzeros)
}


# Calcola la dimensione del file  preso in considerazione
compute_filesize <- function(filename) {
  fileSize <- file.info(filename)$size
  cat("Dimensione del file:", fileSize, "byte\n")
  return(fileSize)
}

process_file <- function(filename, cartella) {
  filename <- gsub(".//", "", filename)
  cat("------------------------------ Elaborazione file", filename, " ------------------------------\n")
  A <- load_matrix_from_file(file.path(cartella, filename))
  if (is_symmetric(A)) {
    cat("La matrice A è simmetrica \n")
  } else {
    cat("La matrice A non è simmetrica \n")
  }
  #DEBUG: print(A)
  b <- create_b_vector(A)
  start_time <- Sys.time()
  resultChol <- cholesky_decomposition(A)
  
  # DEBUG: print(result$memoria_utilizzata_chol)
  
  # DEBUG: print(result$memoria_utilizzata_chol)
  
  resultLin <- solve_linear_system(resultChol$factor, b)
  soluzione <- resultLin[[1]]
  memoria_utilizzata_sistemaLin <- resultLin[[2]]
  
  # DEBUG: print(memoria_utilizzata_sistemaLin)
  memoria_totale<- resultChol$memoria_utilizzata_chol + memoria_utilizzata_sistemaLin
  memoria_totale_mb <- memoria_totale/(1024*1024)
  cat("Memoria totale utilizzata nella risoluzione: ", round(memoria_totale_mb, 2), " MB\n")
  
  end_time <- Sys.time()
  tempo_cholesky <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Stampa il tempo impiegato per la decomposizione di Cholesky
  cat("Tempo di esecuzione per la decomposizione di Cholesky e risoluzione sistema lineare: ", round(tempo_cholesky, 4), " secondi\n")
  
  # DEBUG: print(soluzione)
  # DEBUG: print(memoria_utilizzata_sistemaLin)
  
  # DEBUG: print(soluzione)
  
  errore_relativo <- compute_relative_error(soluzione)
  
  cat("Errore relativo: ", errore_relativo, "\n")
  
  percent_zero <- compute_percent_zeros(A)
  num_nonzeros <- compute_num_nonzeros(A)
  
  fileSize <- compute_filesize(file.path(cartella, filename))
  cat("\n")
  
  return(list(tempo_cholesky = tempo_cholesky, errore_relativo = errore_relativo, percent_zero = percent_zero, 
              num_nonzeros = num_nonzeros, memoria_totale = memoria_totale, fileSize = fileSize, errore_relativo))
}

# # Definiamo cinque vettori per la creazione dei plot
# tempi_totali <- vector()
# memoria_cholesky <- vector()
# nomi_matrici <- character()
# file_size <- numeric()
# errori_relativi <- numeric()


# Crea vettori vuoti per le metriche
tempi_totali <- numeric(0)
memoria_cholesky <- numeric(0)
nomi_matrici <- character(0)
file_size <- numeric(0)
errori_relativi <- numeric(0)

# Cartella contenente i file .mtx
cartella <- "./"

# Esegue il loop dei file contenuti nella cartella
files <- list.files(path = cartella, pattern = "\\.mtx$", full.names = TRUE)

# Ordina i file .mtx in base alla dimensione del file
files <- files[order(file.info(files)$size)]

files <- gsub("\\.//", "", files)
files <- gsub("\"", "", files)

print(files)
# Processa i file .mtx nell'ordine desiderato
for (filename in files) {
  result <- process_file(filename, cartella)
  tempo_cholesky <- result[[1]]
  errore_relativo <- result[[2]]
  memoria_totale <- result[[5]]
  fileSize <- result[[6]]
  
  tempi_totali <- c(tempi_totali, tempo_cholesky)
  
  filename <- gsub("\\.//", "", filename)
  filename <- gsub("\"", "", filename)
  
  nomi_matrici <- c(nomi_matrici, filename)
  memoria_cholesky <- c(memoria_cholesky, memoria_totale)
  file_size <- c(file_size, fileSize)
  errori_relativi <- c(errori_relativi, errore_relativo)
}

# Ordina le metriche in base alla dimensione dei file
sorting_order <- order(file_size)
tempi_totali <- tempi_totali[sorting_order]
memoria_cholesky <- memoria_cholesky[sorting_order]
nomi_matrici <- nomi_matrici[sorting_order]
errori_relativi <- errori_relativi[sorting_order]
file_size <- file_size[sorting_order]

# Ottieni il nome del sistema operativo
operating_system <- Sys.info()["sysname"]

# Imposta il nome del file CSV in base al sistema operativo
if (operating_system == "Linux") {
  filename <- "dati_r_linux.csv"
} else if (operating_system == "Windows") {
  filename <- "dati_r_windows.csv"
} else {
  filename <- "dati_r.csv"  # Nome predefinito nel caso in cui il sistema operativo non sia riconosciuto
}

# Creazione del data frame con le metriche
data <- data.frame(
  MatrixName = nomi_matrici,
  Size = file_size,
  MemoryDiff = memoria_cholesky * 1024 * 1024,  # Converti MB in byte
  Time = tempi_totali,
  Error = errori_relativi
)

# Scrivi il data frame nel file CSV
write.csv(data, file = filename, row.names = FALSE, quote = FALSE)

cat("\nFile CSV creato\n")

