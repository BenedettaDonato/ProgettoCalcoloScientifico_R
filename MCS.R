
# matrix_file <- list.files(path = "/Users/benny/MCS/MatriciMM/", pattern = "\\.mat$", full.names = TRUE)

# Funzione per caricare la matrice dal file .mat e restituire la matrice sparsa
library(Matrix)
library(spam)
library(spam64)
library(tictoc)
library(matrixcalc)
library(pryr)
library(here)



# Funzione per caricare la matrice dal file .mtx e restituire la matrice sparsa
load_matrix_from_file <- function(filepath) {
  A <- readMM(filepath)
  
  # Verifica se la matrice è sparsa
  if (is(A, "CsparseMatrix")) {
    print("La matrice A è sparsa")
  }
  
  return(A)
}

# Definisci il percorso del file .mtx
filepath <- here("MatriciMM", "ex15.mtx")

# Carica e stampa la matrice dal file .mtx
A <- load_matrix_from_file(filepath)

# Controlla se la matrice è simmetrica controllando i valori non nulli della trasposta
is_symmetric <- function(A) {
  A_transpose <- t(A)
  return(all(A@x == A_transpose@x))
}
if (is_symmetric(A)) {
  print("La matrice A è simmetrica")
} else {
  print("La matrice A non è simmetrica")
}


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
  
  # Calcola la fattorizzazione di Cholesky 
  tryCatch({
    factor <- Cholesky(A)
    print("La matrice A è definita positiva")
  }, error = function(err) {
    print("La matrice A non è definita positiva")
  })
  
  # Ottieni la memoria usata dopo la fattorizzazione di Cholesky
  memoria_finale <- mem_used() #usare memory.size() con windows
  
  # Calcola la memoria utilizzata dalla funzione
  memoria_utilizzata_chol <- memoria_finale - memoria_iniziale
  
  return(list(factor = factor, memoria_utilizzata_chol = memoria_utilizzata_chol))
}

# Esegui la fattorizzazione di Cholesky della matrice 'A' e salva il risultato in 'result'
result <- cholesky_decomposition(A)

# Controlla se 'result' contiene la matrice triangolare inferiore 'L'
if (!is.null(result$factor)) {
  # Stampa la matrice triangolare inferiore 'L'
  print(result$factor)
} else {
  print("La matrice A non è definita positiva")
}

print(result$factor) # Matrice triangolare inferiore L (fattorizzazione di Cholesky)
print(result$memoria_utilizzata_chol) # Memoria utilizzata durante la fattorizzazione


# Risolve il sistema lineare Ax=b prendendo in ingresso la matrice 'A' fattorizzata con il metodo di Cholesky e il
# vettore 'b'. Calcola anche la memoria utilizzata durante il processo
solve_linear_system <- function(factor, b) {
  # Ottieni la memoria usata prima della fattorizzazione di Cholesky
  memoria_iniziale <- proc.time()[3]
  
  x <- solve(factor, b)
  
  # DEBUG: print(x)
  
  # Ottieni la memoria usata dopo la fattorizzazione di Cholesky
  memoria_finale <- proc.time()[3]
  
  cat("\n\nSistema Risolto\n")
  
  # Calcola la memoria utilizzata dalla funzione
  memoria_utilizzata_sistemaLin <- memoria_finale - memoria_iniziale
  
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
  cat("------------------------------ Elaborazione file ", filename, " ------------------------------\n")
  A <- load_matrix_from_file(file.path(cartella, filename))
  #DEBUG: print(A)
  b <- create_b_vector(A)
  start_time <- Sys.time()
  result <- cholesky_decomposition(A)
  
  # DEBUG: print(result$memoria_utilizzata_chol)
  
  # DEBUG: print(result$memoria_utilizzata_chol)
  
  result <- solve_linear_system(result$factor, b)
  soluzione <- result[[1]]
  memoria_utilizzata_sistemaLin <- result[[2]]
  
  print(memoria_utilizzata_sistemaLin)
  
  memoria_totale <- result$memoria_utilizzata_chol + memoria_utilizzata_sistemaLin
  
  cat("Memoria totale utilizzata nella risoluzione: ", round(memoria_totale, 2), " MB\n")
  
  end_time <- Sys.time()
  tempo_cholesky <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Stampa il tempo impiegato per la decomposizione di Cholesky
  cat("Tempo di esecuzione per la decomposizione di Cholesky e risoluzione sistema lineare: ", tempo_cholesky, " secondi\n")
  
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

# Definiamo cinque vettori per la creazione dei plot
tempi_totali <- vector()
memoria_cholesky <- vector()
nomi_matrici <- character()
file_size <- numeric()
errori_relativi <- numeric()

# Cartella contenente i file .mtx
cartella <- "./"

# Crea una lista vuota per i file .mtx
file_mat <- list()

# Esegue il loop dei file contenuti nella cartella
files <- list.files(path = cartella, pattern = "\\,mtx$", full.names = TRUE)

# Ordina i file .mtx in base alla dimensione del file
files <- files[order(file.info(files)$size)]

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

# Processa i file .mtx nell'ordine desiderato
for (filename in files) {
  result <- process_file(filename, cartella)
  tempo_cholesky <- result[[1]]
  errore_relativo <- result[[2]]
  memoria_totale <- result[[5]]
  fileSize <- result[[6]]
  
  tempi_totali <- c(tempi_totali, tempo_cholesky)
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

# Crea una lista di etichette per le matrici con dimensioni
matrici_dimensioni <- paste(nomi_matrici, "(", file_size / (1024 * 1024), "MB)", sep="")

# Creazione dei grafici

# Filtra i valori finiti da tempi_totali
valid_tempi_totali <- tempi_totali[is.finite(tempi_totali)]

# Crea il grafico solo se ci sono valori validi
if(length(valid_tempi_totali) > 0) {
  plot(valid_tempi_totali, xlab="Matrici", ylab="Tempo (s)", main="Tempo di Cholesky + Risoluzione per Matrici", 
       xaxt="n", pch=19, col="blue")
  axis(1, at=1:length(nomi_matrici), labels=matrici_dimensioni, las=2)
  abline(h=mean(valid_tempi_totali), col="red")
} else {
  print("Nessun dato valido per creare il grafico.")
}

# Grafico della memoria utilizzata
plot(memoria_cholesky, xlab="Matrici", ylab="Memoria Utilizzata (MB)", main="Memoria utilizzata per Matrici", 
     xaxt="n", pch=19, col="green")
axis(1, at=1:length(nomi_matrici), labels=matrici_dimensioni, las=2)
abline(h=mean(memoria_cholesky), col="red")

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
write.csv(data, file = filename, row.names = FALSE)

cat("\nFile CSV creato\n")

