#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iostream"
#include "random"
#include <pthread.h>

#include "chrono"
#include "time.h"
#include "cassert"

using namespace std;
using namespace std::chrono;

int n; // Size of the matrix
int n2; // Size of the matrix squared
double *U, *L, *a, *original_a;
int *pi;
int global_seed;
double max_val;
int max_idx;
int num_threads;
pthread_t *threads;

typedef struct {
    int start; // Start index of the iteration range
    int end;   // End index of the iteration range
	unsigned int rand_seed;
	int id;
	int k;
	int row_offset;
    pthread_mutex_t *mutex;
} ThreadArgs;

ThreadArgs *thread_data;

// Function to initialize matrices and arrays
void *initialize_matrices(void* arg) {
	ThreadArgs* thread_data = (ThreadArgs*)arg;
	int i, j, k;
	// generate random matrix
	for (k = thread_data->start; k < thread_data->end; k++) {
		i = k / n;
		j = k % n;
		a[k] = rand_r(&thread_data->rand_seed)/(double)RAND_MAX;
		original_a[k] = a[k];
		if (i > j)
			U[k] = 0; // 0's below the diagonal
		else if (i == j)
			L[k] = 1; // 1's on the diagonal
		else
			L[k] = 0; // 0's above the diagonal
		if(j==0)
			pi[i] = i;	// identity permutation
	}
	return NULL;
}

void *find_max(void *arg) {
    ThreadArgs *data = (ThreadArgs *)arg;

    double local_max = 0.0;
    int local_max_idx = -1;

    for (int i = data->start; i < data->end; i++) {
        if (local_max < fabs(a[i * n + data->k])) {
            local_max = fabs(a[i * n + data->k]);
            local_max_idx = i;
        }
    }

    pthread_mutex_lock(data->mutex);
    if (local_max > max_val) {
        max_val = local_max;
        max_idx = local_max_idx;
    }
    pthread_mutex_unlock(data->mutex);

    pthread_exit(NULL);
}

void pthread_find_max(int k){
	max_idx = -1;
	max_val = 0.0;

	pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	int n_minus_k = n - k;
	for (int i = 0; i < num_threads; i++) {
		thread_data[i].mutex = &mutex;
		thread_data[i].k = k;
		thread_data[i].start = (i * n_minus_k) / num_threads + k;
		thread_data[i].end = ((i + 1) * n_minus_k) / num_threads + k;

		pthread_create(&threads[i], NULL, find_max, (void *)&thread_data[i]);
	}

	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);

	if (max_val == 0) {// comparisons among float is bad, replace with |max| < eps ?
		cerr << "Singular Matrix\n";
		return;
	}

	swap(pi[k], pi[max_idx]);
}

void *swap_rows(void *arg) {
    ThreadArgs *data = (ThreadArgs *)arg;
	int row_offset = data->k * n;
    for (int j = data->start; j < data->end; j++) {
        swap(a[row_offset + j], a[max_idx * n + j]);
        if (j < data->k) {
            swap(L[row_offset + j], L[max_idx*n + j]);
        }
    }
    pthread_exit(NULL);
}

void pthread_swap_rows(int k){
	for (int i = 0; i < num_threads; i++) {
		thread_data[i].k = k;
		thread_data[i].start = (i * n) / num_threads;
		thread_data[i].end = ((i + 1) * n) / num_threads;

		pthread_create(&threads[i], NULL, swap_rows, (void *)&thread_data[i]);
	}

	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);
}

void *update_matrices(void *arg){
	ThreadArgs *data = (ThreadArgs *)arg;
	int row_offset = data->row_offset;
	for (int i = data->start; i < data->end; i++) {
		L[i*n + data->k] = a[i*n + data->k]/U[row_offset + data->k];
		U[row_offset + i] = a[row_offset + i];
		for (int j = data->k + 1; j < n; j++)
			a[i*n + j] -= L[i*n + data->k] * a[row_offset + j];
	}
	pthread_exit(NULL);

}

void pthread_update_matrices(int k){
	int row_offset = k * n;
	U[row_offset + k] = a[row_offset + k];
	int n_minus_k_minus_1 = n - k - 1;
	for (int i = 0; i < num_threads; i++) {
		thread_data[i].k = k;
		thread_data[i].row_offset = row_offset;
		thread_data[i].start = (i * n_minus_k_minus_1) / num_threads + k+1;
		thread_data[i].end = ((i + 1) * n_minus_k_minus_1) / num_threads + k+1;

		pthread_create(&threads[i], NULL, update_matrices, (void *)&thread_data[i]);
	}

	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);
}

void decompose () {

	for (int k = 0; k < n; k++) {
		
		pthread_find_max(k);

		pthread_swap_rows(k);
		
		pthread_update_matrices(k);
		
	}
}

void print(){
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cout << original_a[i*n + j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	for (i = 0; i < n; i++){
		cout << pi[i] << " ";
	}
	cout << "\n\n";

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cout << L[i*n + j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cout << U[i*n + j] << " ";
		}
		cout << "\n";
	}
}

int main (int argc, char* argv[]) {

	n = atoi(argv[1]);
	n2 = n * n;
	num_threads = atoi(argv[2]);
	thread_data = (ThreadArgs *)malloc(num_threads * sizeof(ThreadArgs));
	threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

	global_seed = time(NULL);

	// Initialize random number generator for each thread
    for (int i = 0; i < num_threads; i++) { // Adjust the number of threads accordingly
		thread_data[i].id = i;
    	thread_data[i].rand_seed = global_seed + i;
    }

	U = (double *)malloc(n2 * sizeof(double));
    L = (double *)malloc(n2 * sizeof(double));
	a = (double *)malloc(n2 * sizeof(double));
	original_a = (double *)malloc(n2 * sizeof(double));
    pi = (int *)malloc(n * sizeof(int));
	
	 // Create threads
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].start = (i * n2) / num_threads;
        thread_data[i].end = ((i + 1) * n2) / num_threads;
        pthread_create(&threads[i], NULL, initialize_matrices, (void *)&thread_data[i]);
    }

	// Join threads
    for (int i = 0; i < num_threads; i++)
        pthread_join(threads[i], NULL);

    auto start = high_resolution_clock::now();
    decompose();
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<nanoseconds>(stop - start);
    cerr << duration.count() << endl;




    print();


	return 0;
}