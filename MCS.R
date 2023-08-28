invisible({
  # Funzione per caricare la matrice dal file .mat e restituire la matrice sparsa
  library(Matrix)
  #library(spam)
  #library(spam64)
  library(tictoc)
  library(matrixcalc)
  library(pryr)
  library(here)
  #library(SparseM)
  library(R.matlab)
  
  
  # Carica le funzioni dal file delle funzioni
  source("Utils.R")
  
  # Crea vettori vuoti per le metriche
  tempi_totali <- numeric(0)
  memoria_cholesky <- numeric(0)
  nomi_matrici <- character(0)
  file_size <- numeric(0)
  errori_relativi <- numeric(0)
  
  # Cartella contenente i file .mat 
  cartella <- here()
  
  # Esegue il loop dei file contenuti nella cartella
  files <- list.files(path = cartella, pattern = "\\.mat$", full.names = FALSE)
  
  # Ordina i file .mat in base alla dimensione del file
  files <- files[order(file.info(files)$size)]
  
  # Si pulisce il nome del file .mat 
  files <- gsub("\\.//", "", files)
  files <- gsub("\"", "", files)
  
  # Processa i file .mat nell'ordine desiderato
  for (filename in files) {
    result <- process_file(filename, cartella)
    tempo_cholesky <- result[[1]]
    errore_relativo <- result[[2]]
    memoria_totale <- result[[5]]
    fileSize <- result[[6]]
    
    tempi_totali <- c(tempi_totali, tempo_cholesky)
    
    # Eliminiamo le virgolette e il numero prima dei file .mat
    filename <- gsub("\\.//", "", filename)
    filename <- gsub("\"", "", filename)
    filename <- gsub("^[0-9 ]+", "", filename)
    
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
    MemoryDiff = memoria_cholesky,
    Time = tempi_totali,
    Error = errori_relativi
  )
  
  # Scrivi il data frame nel file CSV
  write.csv(data, file = filename, row.names = FALSE, quote = FALSE)
  
  cat("\nFile CSV creato\n")
  
  
})
  


