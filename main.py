import numpy as np

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):

    temp_seq1 = seq1
    # Inicialización de la matriz de puntuaciones
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int)
    max_score = 0
    max_positions = []

    # Llenado de la matriz
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch

            score_matrix[i][j] = max(
                score_matrix[i - 1][j - 1] + score,  # Diagonal
                score_matrix[i - 1][j] + gap,        # Arriba
                score_matrix[i][j - 1] + gap,        # Izquierda
                0                                    # Cero
            )

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_positions = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_positions.append((i, j))

    # Trazado hacia atrás para encontrar la subsecuencia de puntaje máximo
    alignments = []
    for pos in max_positions:
        i, j = pos
        aligned_seq = []
        end_i, end_j = i, j
        while score_matrix[i][j] != 0:
            aligned_seq.append(seq1[i - 1])
            if score_matrix[i][j] == score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
                if score_matrix[i][j] > score_matrix[i - 1][j - 1]:
                    i -= 1
                    j -= 1
                else:
                    aligned_seq.pop()
                    break
            else:
                break
        aligned_seq.reverse()
        alignments.append((''.join(seq1), ''.join(seq2), ''.join(aligned_seq), end_i + 1, end_j + 1, i + 1, j + 1))

    return score_matrix, max_score, alignments

def save_to_file(filename, score_matrix, max_score, alignments):
    with open(filename, 'w') as f:
        f.write("Matriz de scores:\n")
        for row in score_matrix:
            f.write(' '.join(map(str, row)) + '\n')
        
        f.write(f"\nPuntaje maximo: {max_score}\n")
        
        for idx, (seq1, seq2, subseq, end_pos1, end_pos2, start_pos1, start_pos2) in enumerate(alignments, start=1):
            f.write(f"\nSubsecuencia {idx} del puntaje maximo:\n")
            f.write(f"Subsecuencia: {subseq}\n")
            f.write(f"Secuencia 1: {seq1}\n")
            f.write(f"Secuencia 2: {seq2}\n")
            f.write(f"Posicion en Secuencia 1: {start_pos1} a {end_pos1}\n")
            f.write(f"Posicion en Secuencia 2: {start_pos2} a {end_pos2}\n")

# Ejemplo de uso
seq1 = "GTACACGT"
seq2 = "GTTCACG"

seq3 = "ACGTACTA"
seq4 = "GTCGTAG"

seq3 = "acgatagcagatagcgcatagcgactagcgactgcagctacgcagcatagcagcagcagaacgatagcagatagcgcatagcgactagcgactgcagctacgcagcatagcagcagcaga"
seq4 = "tgagctagagatagctacgacgcatcagcgatagcagctaggcagctgcagcgactagcatgagctagagatagctacgacgcatcagcgatagcagctaggcagctgcagcgactagca"

# seq3 = "AACGTTTACA"
# seq4 = "AGCACACA"

# score_matrix, max_score, alignments = smith_waterman(seq1, seq2)
score_matrix, max_score, alignments = smith_waterman(seq3, seq4)
save_to_file("smith_waterman_output.txt", score_matrix, max_score, alignments)
