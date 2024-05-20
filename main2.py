import numpy as np

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
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

            diagonal_score = score_matrix[i - 1][j - 1] + score
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diagonal_score, up_score, left_score, 0)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_positions = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_positions.append((i, j))

    # Trazado hacia atrás para encontrar la subsecuencia de puntaje máximo
    alignments = []
    for pos in max_positions:
        i, j = pos
        aligned_seq1 = []
        aligned_seq2 = []
        while score_matrix[i][j] != 0:
            if score_matrix[i][j] == score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_matrix[i - 1][j] + gap:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1
        
        # Reverse aligned sequences to get correct order
        aligned_seq1.reverse()
        aligned_seq2.reverse()

        # Filter out consecutive identical characters
        filtered_seq1 = ''
        filtered_seq2 = ''
        for idx in range(len(aligned_seq1)):
            if aligned_seq1[idx] == aligned_seq2[idx]:
                filtered_seq1 += aligned_seq1[idx]
                filtered_seq2 += aligned_seq2[idx]

        alignments.append((filtered_seq1, filtered_seq2, pos[0], pos[1]))

    return score_matrix, max_score, alignments

def save_to_file(filename, score_matrix, max_score, alignments):
    with open(filename, 'w') as f:
        f.write("Matriz de scores:\n")
        for row in score_matrix:
            f.write(' '.join(map(str, row)) + '\n')
        
        f.write(f"\nPuntaje máximo: {max_score}\n")
        
        for idx, (seq1, seq2, end_pos1, end_pos2) in enumerate(alignments, start=1):
            f.write(f"\nSubsecuencia {idx} del puntaje máximo:\n")
            f.write(f"Secuencia 1: {seq1}\n")
            f.write(f"Secuencia 2: {seq2}\n")
            f.write(f"Posición en Secuencia 1: {end_pos1 - len(seq1) + 1} a {end_pos1}\n")
            f.write(f"Posición en Secuencia 2: {end_pos2 - len(seq2) + 1} a {end_pos2}\n")

# Ejemplo de uso
seq1 = "ACCGT"
seq2 = "ACG"

seq3 = "acgatagcagatagcgcatagcgactagcgactgcagctacgcagcatagcagcagcagaacgatagcagatagcgcatagcgactagcgactgcagctacgcagcatagcagcagcaga"
seq4 = "tgagctagagatagctacgacgcatcagcgatagcagctaggcagctgcagcgactagcatgagctagagatagctacgacgcatcagcgatagcagctaggcagctgcagcgactagca"

# score_matrix, max_score, alignments = smith_waterman(seq1, seq2)
score_matrix, max_score, alignments = smith_waterman(seq3, seq4)
save_to_file("smith_waterman_output.txt", score_matrix, max_score, alignments)
