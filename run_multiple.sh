#!/bin/bash

# Parametri iniziali
start_ind=0
end_ind=200
increment=200  # Incremento per start_ind e end_ind
output_prefix="out_cal"

# Numero di iterazioni
num_iterations=100

# Ciclo
for ((i=1; i<=num_iterations; i++))
do
  # Modifica i parametri nel file di configurazione
  sed -i "s/^start_ind=.*/start_ind=${start_ind}/" config.txt
  sed -i "s/^end_ind=.*/end_ind=${end_ind}/" config.txt
  sed -i "s/^output_file=.*/output_file=${output_prefix}_${i}.txt/" config.txt

  # Esegui il programma
  ./runfit config.txt

  # Attendi che il programma finisca prima di continuare
  wait
  
  # Aggiorna i valori per la prossima iterazione
  start_ind=$((start_ind + increment))
  end_ind=$((end_ind + increment))
  

done
