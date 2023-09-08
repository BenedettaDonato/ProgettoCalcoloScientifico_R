# ProgettoCalcoloScientifico_R
> Progetto R per la risoluzione con il metodo di Cholesky di sistemi lineari con matrici sparse e definite positive di grandi dimensioni.

## Introduzione

Questo Repository contiene il programma scritto in R necessario per completare il progetto di Metodi del Calcolo Scientifico. 
All'interno del progetto vengono utiliizzate le librerie esterne, per la risoluzione di un sistema lineare associato ad una matrice sparsa, simmetrica e definita positiva utilizzando la decomposizione di Cholesky.

## Documentazione utilizzata
Link utile per visionare la documentazione delle librerie utilizzate durante l'implementazione del progetto:
- Documentazione pryr: [https://docs.scipy.org/doc/scipy/](https://cran.r-project.org/web/packages/pryr/)
- Documentazione here: [https://scikit-sparse.readthedocs.io/en/latest/](https://cran.r-project.org/web/packages/here/)
- Documentazione R.matlab: [https://scikit-sparse.readthedocs.io/en/latest/](https://cran.r-project.org/web/packages/R.matlab/index.html)

## Note Generali
Il file *MCS.R* contiene il codice per eseguire la decomposizione di Cholesky e la risoluzione 
del sistema lineare *Ax=b* per gran parte delle matrici richieste dalla consegna. 

## Matrici testate
All'interno del programma sono state testate tutte le matrici richieste dalla specifica del progetto. Tuttavia, non tutte vengono supportate dalle librerie da noi utilizzate in quanto quelle di dimensioni maggiori o con una percentuale di elementi non nulli alta, esauriscono velocemente le risorse del sistema.
Qui di seguito è mostrata una tabella che elenca tutte le matrici da noi testate affiancate da un indicatore indicante l’esito (positivo o negativo) del processo di fattorizzazione e risoluzione del sistema lineare: 

| Matrici | Stato | Dimensione (in MB) |
|------------------------|------------------------|-------------------|
| [ex15.mat](https://sparse.tamu.edu/FIDAP/ex15)   | ✅ | 0,555        |
| [shallow_water1.mat](https://sparse.tamu.edu/MaxPlanck/shallow_water1) | ✅ |   2,263      |
| [apache2.mat](https://sparse.tamu.edu/GHS_psdef/apache2) | ✅ | 8,302      |
| [parabolic_fem.mat](https://sparse.tamu.edu/Wissgott/parabolic_fem)  | ✅ | 13,116       |
| [G3_circuit.mat](https://sparse.tamu.edu/AMD/G3_circuit) | ✅ | 13,833      |
| [cfd1.mat](https://sparse.tamu.edu/Rothberg/cfd1) | ✅ | 14,164      |
| [cfd2.mat](https://sparse.tamu.edu/Rothberg/cfd2) | ✅ | 23,192       |
| [StocF - 1465.mat](https://sparse.tamu.edu/Janna/StocF-1465)  | ❌ Out of memory | 178,368        |
| [Flan_1565.mat](https://sparse.tamu.edu/Janna/Flan_1565)  | ❌ Out of memory | 292,858  |

## Descrizione del programma
Il programma esegue diverse operazioni su una serie di file che contengono matrici sparse e definite positive, con lo scopo di calcolare il tempo di esecuzione, la memoria utilizzata e l'errore relativo, tra la soluzione calcolata e quella esatta, durante la decomposizione di Cholesky e la risoluzione di un sistema lineare *Ax=b*. 



## Librerie utilizzate
La libreria pryr, disponibile in R, è uno strumento per gli sviluppatori che offre un insieme di funzioni e strumenti per il monitoraggio e l'analisi della memoria durante l’esecuzione di calcoli complessi, consentendo agli utenti di diagnosticarla e ottimizzarla nei loro script.
È particolarmente adatta quando si lavora con grandi insiemi di dati in quanto aiuta a individuare potenziali problemi di memoria e fornisce informazioni utili a migliorare l'efficienza del codice. Inoltre, questa libreria svolge un ruolo importante nel garantire che i programmi siano ottimizzati per l'esecuzione su sistemi con risorse limitate, migliorando così la riproducibilità e l'efficienza delle analisi scientifiche condotte.

La libreria here, importata nello script R, è uno strumento essenziale per semplificare la gestione delle directory all'interno dei progetti. 
here può essere descritta come uno strumento fondamentale per migliorare la riproducibilità e la portabilità delle analisi condotte; di fatto facilita la gestione dei percorsi all'interno dei progetti, semplifica la condivisione del codice e migliora la trasferibilità delle analisi tra sistemi operativi. 
Questo contribuisce notevolmente ad aumentare la robustezza e l'affidabilità delle ricerche scientifiche condotte, assicurando che gli altri utilizzatori possano facilmente riprodurre i risultati ottenuti.

La libreria R.matlab è uno strumento fondamentale per l'interazione e lo scambio di dati tra R e il formato file MATLAB (.mat). 
R.matlab può essere descritta come una libreria essenziale per l'integrazione di dati provenienti da software e ambienti di calcolo diversi offrendo “un ponte” tra R e MATLAB e consentendo agli utenti di sfruttare al meglio le potenzialità di entrambi i linguaggi. 

Oltre alle librerie illustrate precedentemente, ne sono state utilizzate delle altre, disponibili nativamente in R, per l'analisi e la manipolazione dei dati; così facendo, è stato possibile sfruttare le potenzialità del linguaggio R semplificando la gestione dei dati provenienti da diverse fonti, condurre analisi avanzate e risultati confrontabili.
