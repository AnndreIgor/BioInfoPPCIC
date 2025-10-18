from Bio import SeqIO
import os

def analisar_fasta(caminho_fasta):
    """
    Analisa um arquivo FASTA e retorna estatísticas básicas sobre as sequências.

    Parâmetros:
    -----------
    caminho_fasta : str
        Caminho do arquivo FASTA a ser analisado.

    Retorna:
    --------
    dict
        Dicionário contendo:
        - 'quantidade': número total de sequências
        - 'maior': comprimento da maior sequência
        - 'menor': comprimento da menor sequência
        - 'media': comprimento médio das sequências
    """
    
    comprimentos = [len(record.seq) for record in SeqIO.parse(caminho_fasta, "fasta")]
    
    if not comprimentos:
        return {
            "quantidade": 0,
            "maior": 0,
            "menor": 0,
            "media": 0.0
        }

    return {
        "arquivo": os.path.basename(caminho_fasta),
        "quantidade": len(comprimentos),
        "maior": max(comprimentos),
        "menor": min(comprimentos),
        "media": round(sum(comprimentos) / len(comprimentos), 2),
        "size_bytes": os.path.getsize(caminho_fasta)
    }

# Exemplo de uso:
if __name__ == "__main__":
    caminho = r"data/full_dataset_plasmodium/PLASMODIUM0.fasta"  # Substitua pelo caminho do seu arquivo
    resultado = analisar_fasta(caminho)
    print(resultado)