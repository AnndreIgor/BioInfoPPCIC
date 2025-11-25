#!/bin/bash
#$ -cwd                    # Executar no diretório atual
#$ -o saida.log            # Arquivo de saída padrão
#$ -e erro.log             # Arquivo de erro
#$ -V                      # Exportar variáveis de ambiente

echo "Iniciando job em $(hostname) com $NSLOTS núcleos"
echo "Diretório de execução: $(pwd)"
date

# Aqui entra o comando Python
cd code
../.venv/bin/python main.py

date
echo "Job finalizado."
