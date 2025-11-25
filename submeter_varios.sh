#!/bin/bash
# Script para submeter várias vezes o mesmo job com núcleos aleatórios

NUM_EXEC=100   # Quantas vezes rodar o job (ajuste aqui)
JOB="job_python.sh"  # Script do job que será submetido

for ((i=1; i<=NUM_EXEC; i++)); do
    # Gera número aleatório entre 1 e 24
    NCORES=$(( ( RANDOM % 24 ) + 1 ))

    echo "[$i/$NUM_EXEC] Submetendo com $NCORES núcleos..."
    qsub -N "job_${NCORES}cores" -pe smp $NCORES $JOB
done
