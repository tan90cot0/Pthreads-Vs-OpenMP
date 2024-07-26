import subprocess
import numpy as np
from scipy.linalg import lu

def run_LU_serial(n):
    compile_process = subprocess.Popen(["make", "serial"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
        process = subprocess.Popen(["./LU_serial.o", str(n), str(1)], stdout=subprocess.PIPE, universal_newlines=True)
        output = process.communicate()[0]
        return output.split("\n")
    else:
        print('Compilation failed')

def run_LU_pthreads(n, t):
    compile_process = subprocess.Popen(["make", "pthreads"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
        process = subprocess.Popen(["./LU_pthreads.o", str(n), str(t)], stdout=subprocess.PIPE, universal_newlines=True)
        output = process.communicate()[0]
        return output.split("\n")
    else:
        print('Compilation failed')

def run_LU_openmp(n, t):
    # Compile the program using 'make openmp'
    compile_process = subprocess.Popen(["make", "openmp"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
    # Run the compiled program with specified arguments
        process = subprocess.Popen(["./LU_openmp.o", str(n), str(t)], stdout=subprocess.PIPE, universal_newlines=True)
        output, _ = process.communicate()
        return output.split("\n")

    else:
        print('Compilation failed')
def parse_matrix(rows):
    return np.array([list(map(float, row.split())) for row in rows])

def parse_permutation_matrix(permutation_string):
    permutation = list(map(int, permutation_string.split()))
    n = len(permutation)
    identity_matrix = np.eye(n)
    return identity_matrix[permutation]

def check_PALU(A, P, L, U):
    PA = np.matmul(P,A)
    LU = np.matmul(L,U)
    # print(PA - LU)
    return np.allclose(PA, LU, atol=1e-3)

def general_checker(n, output):
    output = [s for s in output if s]
    A = parse_matrix(output[:n])
    P = parse_permutation_matrix(output[n])
    L = parse_matrix(output[n+1: 2*n + 1])
    U = parse_matrix(output[2*n + 1: 3*n + 1])

    P_c, L_c, U_c = lu(A)
    P_c = np.linalg.inv(P_c)

    tolerance = 1e-2

    if np.allclose(U, U_c, atol=tolerance):
        print("U is correct")

    if np.allclose(L, L_c, atol=tolerance):
        print("L is correct")

    if np.allclose(P, P_c, atol=tolerance):
        print("P is correct")

if __name__ == "__main__":
    n = 2000  # Change this to the desired size
    t = 6  # Change this to the desired number of threads
    output_pthreads = run_LU_pthreads(n, t)
    general_checker(n, output_pthreads)

    # for t in [3]:
    #     output_openmp = run_LU_openmp(n, t)
    #     general_checker(n, output_openmp)


