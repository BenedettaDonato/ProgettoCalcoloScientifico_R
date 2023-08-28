# File R contenente le funzioni principali per elaborare le matrici in formato .mat

# Funzione per caricare la matrice dal file .mat e restituire la matrice sparsa
load_matrix_from_file <- function(filepath) {
  # A <- readMM(filepath)
  #DEBUG: cat(filepath)
  mat_data <- readMat(filepath)
  
  problem_elements <- mat_data$Problem
  
  
  A <- NULL
  
  for (element in problem_elements) {
    if (inherits(element, "dgCMatrix")) {
      A <- element
      break  # Esci dal loop una volta trovato l'oggetto
    }
  }
  
  # Verifica se la matrice è sparsa
  if (is.sparseMatrix(A)) {
    cat("La matrice A è sparsa \n")
  }
  return(A)
}

is.sparseMatrix <- function(x) is(x, 'sparseMatrix')


# Definisci il percorso del file .mat
# filepath <- here("MatriciMM", "ex15.mat")

# Carica e stampa la matrice dal file .mat
# A <- load_matrix_from_file(filepath)

# Controlla se la matrice è simmetrica controllando i valori non nulli della trasposta

is_symmetric <- function(A) {
  A_transpose <- t(A)
  return(all(A@x == A_transpose@x))
}


# print(is_symmetric(matrix_A))

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
  memoria_iniziale <- mem_used()
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

# b <- create_b_vector(matrix_A)
# factor_prova <- cholesky_decomposition(matrix_A)



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
  #DEBUG: cat(filename, "FILENAME")
  A <- load_matrix_from_file(file.path(cartella, filename))
  if (is_symmetric(A)) {
    cat("La matrice A è simmetrica \n")
  } else {
    cat("La matrice A non è simmetrica \n")
  }
  # DEBUG: print(A)
  b <- create_b_vector(A)
  # DEBUG: start_time <- Sys.time()
  # DEBUG: cat("Tempo iniziale: ", start_time, "\n")
  start_time_chol <- Sys.time()
  resultChol <- cholesky_decomposition(A)
  end_time_chol <- Sys.time()
  differenza_chol_time <- as.numeric(difftime(end_time_chol, start_time_chol, units = "secs"))
  # DEBUG: print(result$memoria_utilizzata_chol)
  
  # DEBUG: print(result$memoria_utilizzata_chol)
  start_time_lin <- Sys.time()
  resultLin <- solve_linear_system(resultChol$factor, b)
  end_time_lin<- Sys.time()
  differenza_lin_time <- as.numeric(difftime(end_time_lin, start_time_lin, units = "secs"))
  soluzione <- resultLin[[1]]
  memoria_utilizzata_sistemaLin <- resultLin[[2]]
  
  # DEBUG: print(memoria_utilizzata_sistemaLin)
  memoria_totale<- resultChol$memoria_utilizzata_chol + memoria_utilizzata_sistemaLin
  memoria_totale_mb <- memoria_totale/(1024*1024)
  cat("Memoria totale utilizzata nella risoluzione: ", round(memoria_totale_mb, 2), " MB\n")
  
  # DEBUG: end_time <- Sys.time()
  #DEBUG: cat("Tempo finale: ", end_time, "\n")
  tempo_cholesky <- differenza_lin_time + differenza_chol_time
  
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

