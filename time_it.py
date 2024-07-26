import subprocess
import numpy as np
from scipy.linalg import lu
import csv
import itertools
import os
import random


def run_LU_serial(n):
    compile_process = subprocess.Popen(["make", "serial"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
        process = subprocess.Popen(["./LU_serial.o", str(n), str(1)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        _ , timing = process.communicate()
        print(f'Got time for {n} as {timing}')
        return float(timing)
    else:
        print('Compilation failed')

def run_LU_openmp(n, t):
    # Compile the program using 'make openmp'
    compile_process = subprocess.Popen(["make", "openmp"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
        # Run the compiled program with specified arguments
        process = subprocess.Popen(["./LU_openmp.o", str(n), str(t)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        _, timing = process.communicate()
        print(f'Got time for {n},{t} as {timing}')
        return float(timing)

    else:
        print('Compilation failed')

def run_LU_pthreads(n, t):
    # Compile the program using 'make openmp'
    compile_process = subprocess.Popen(["make", "pthreads"], stdout=subprocess.PIPE, universal_newlines=True)
    compile_output, _ = compile_process.communicate()
    print(compile_output)

    if compile_process.returncode == 0:
        # Run the compiled program with specified arguments
        process = subprocess.Popen(["./LU_pthreads.o", str(n), str(t)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        _, timing = process.communicate()
        print(f'Got time for {n},{t} as {timing}')
        return float(timing)

    else:
        print('Compilation failed')


if __name__ == "__main__":
    ns = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    ts = [1, 2, 4, 8, 16, 32, 64]
    combinations = list(itertools.product(ns, ts))
    random.shuffle(combinations)

    openmp_filename = 'output_openmp.csv'
    pthreads_filename = 'output_pthreads.csv'
    serial_filename = 'output_serial.csv'


    existing_openmp_combinations = []
    if os.path.isfile(openmp_filename):
        with open(openmp_filename, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                n = int(row[0])
                t = int(row[1])
                existing_openmp_combinations.append((n, t))

    existing_pthreads_combinations = []
    if os.path.isfile(pthreads_filename):
        with open(pthreads_filename, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                n = int(row[0])
                t = int(row[1])
                existing_pthreads_combinations.append((n, t))

    existing_serial_combinations = []
    if os.path.isfile(serial_filename):
        with open(serial_filename, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                n = int(row[0])
                t = int(row[1])
                existing_serial_combinations.append((n, t))


    for n, t in combinations:

            if (n, t) not in existing_openmp_combinations:
                try:
                    print(f'Running for {n} and {t} with OpenMP')
                    timing_run = run_LU_openmp(n, t)
                    print(f'Got time for {n},{t} as {timing_run}')
                    with open(openmp_filename, 'a') as openmp_file:
                        openmp_writer = csv.writer(openmp_file)
                        openmp_writer.writerow([n, t, timing_run])
                except Exception as e:
                    print(f'Failed for {n},{t} with OpenMP {e}')
                except:
                    print(f'Failed for {n},{t} with OpenMP')

            if (n, t) not in existing_pthreads_combinations:
                try:
                    print(f'Running for {n} and {t} with Pthreads')
                    timing_run = run_LU_pthreads(n, t)
                    print(f'Got time for {n},{t} as {timing_run}')

                    with open(pthreads_filename, 'a') as pthreads_file:
                        pthreads_writer = csv.writer(pthreads_file)
                        pthreads_writer.writerow([n, t, timing_run])
                except Exception as e:
                    print(f'Failed for {n},{t} with Pthreads {e}')
                except:
                    print(f'Failed for {n},{t} with Pthreads')

            if t == 1:
                if (n, t) not in existing_serial_combinations:
                    try:
                        print(f'Running for {n} for serial')
                        timing_run = run_LU_serial(n)
                        print(f'Got time for {n} as {timing_run}')

                        with open(serial_filename, 'a') as serial_file:
                            serial_writer = csv.writer(serial_file)
                            serial_writer.writerow([n, t, timing_run])
                    except Exception as e:
                        print(f'Failed for {n} with serial {e}')
                    except:
                        print(f'Failed for {n} for serial')
