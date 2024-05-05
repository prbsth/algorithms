import sys
import time
import random
import numpy as np

def generate_matrix(dim):
    return np.random.randint(0, 10, size=(dim, dim))

def save_matrices_to_file(A, B, filename):
    with open(filename, 'w') as fout:
        np.savetxt(fout, A, fmt='%d')
        fout.write('\n')
        np.savetxt(fout, B, fmt='%d')

def strassen(A, B, dim, threshold):
    if dim == 1:
        return A * B
    
    if dim <= threshold:
        return standard(A, B)
    
    # add 0 padding if dim is odd
    padded = False
    if dim % 2 == 1:
        padded = True
        dim += 1
        A = np.pad(A, ((0, 1), (0, 1)), mode='constant')
        B = np.pad(B, ((0, 1), (0, 1)), mode='constant')

    mid = dim // 2
    a, b, c, d = A[:mid, :mid], A[:mid, mid:], A[mid:, :mid], A[mid:, mid:]
    e, f, g, h = B[:mid, :mid], B[:mid, mid:], B[mid:, :mid], B[mid:, mid:]

    P1 = strassen(a, f - h, mid, threshold)
    P2 = strassen(a + b, h, mid, threshold)
    P3 = strassen(c + d, e, mid, threshold)
    P4 = strassen(d, g - e, mid, threshold)
    P5 = strassen(a + d, e + h, mid, threshold)
    P6 = strassen(b - d, g + h, mid, threshold)
    P7 = strassen(a - c, e + f, mid, threshold)

    AEBG = P4 + P2 - P5 + P6
    AFBH = P1 + P2
    CEDG = P3 + P4
    CFDH = P1 - P3 + P5 + P7

    result = np.zeros((dim, dim), dtype=int)
    result[:mid, :mid] = AEBG
    result[:mid, mid:] = AFBH
    result[mid:, :mid] = CEDG
    result[mid:, mid:] = CFDH

    if padded:
        result = result[:-1, :-1]

    return result

def standard(A, B):
    return A.dot(B)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python strassen.py <dimension>")
        sys.exit(1)

    max_dimension = int(sys.argv[1])
    dimensions = [2**i for i in range(2, int(np.log2(max_dimension)) + 1)]  # Adjust to go up to max_dimension
    dimensions.append(max_dimension)  # Ensure max_dimension is included if it's not a power of 2

    thresholds = [2**i for i in range(1, 8)]  # Threshold values from 2 to 128

    print("Dimension\tThreshold\tStrassen\tStandard")
    for dimension in dimensions:
        A = generate_matrix(dimension)  # Regenerate A for the current dimension
        B = generate_matrix(dimension)  # Regenerate B for the current dimension

        for threshold in thresholds:
            if threshold > dimension:
                break

            start_time = time.time()
            strassen_result = strassen(A, B, dimension, threshold)
            strassen_time = time.time() - start_time

            start_time = time.time()
            standard_result = standard(A, B)
            standard_time = time.time() - start_time

            print(f"{dimension}\t\t{threshold}\t\t{strassen_time:.6f}\t{standard_time:.6f}")

