#!/bin/bash
# Script: rodar_python_n_vezes.sh

N=3  # número de vezes que o código vai rodar
    
cd code

for i in $(seq 1 $N)
do
    echo "Execução número $i"
    python3 main.py
done

echo "Finalizado após $N execuções."