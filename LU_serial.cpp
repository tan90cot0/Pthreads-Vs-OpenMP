#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iostream"
#include "random"
#include <time.h>

#include "cassert"
#include "chrono"

using namespace std;
using namespace std::chrono;

void decompose (double* a, int* pi, double* L, double* U, int n) {
	int i, j, k;
	// initialize L and U
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i > j)
				U[i*n + j] = 0; // 0's below the diagonal
			else if (i == j)
				L[i*n + j] = 1; // 1's on the diagonal
			else
				L[i*n + j] = 0; // 0's above the diagonal
		}
		pi[i] = i;	// identity permutation
	}

	for (k = 0; k < n; k++) {
		double max = 0;
		int max_idx = 0;
		for (i = k; i < n; i++) {
			if (max < fabs(a[i*n + k])) {
				max = fabs(a[i*n + k]);
				max_idx = i;
			}
		}

        // cerr << "Max: " << max << endl;


        if (max == 0) {// comparisons among float is bad, replace with |max| < eps ?
			// cerr << "Singular Matrix\n";
			return;
		}

        // cerr << "Printing the first A"<< endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                // cerr << a[i * n + j] << " ";
//            }
//            // cerr << endl;
//        }

		swap(pi[k], pi[max_idx]);

		for (j = 0; j < n; j++) {	// hmm there must be an STL for this
			swap(a[k*n + j], a[max_idx*n + j]);
			if (j < k)
				swap(L[k*n + j], L[max_idx*n + j]);
		}

        // cerr << "Printing the second A"<< endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                // cerr << a[i * n + j] << " ";
//            }
//            // cerr << endl;
//        }

		U[k*n + k] = a[k*n + k];

		for (i = k + 1; i < n; i++) {
			L[i*n + k] = a[i*n + k]/U[k*n + k];
			U[k*n + i] = a[k*n + i];
		}

		for (i = k + 1; i < n; i++) {
			for (j = k + 1; j < n; j++) {
				a[i*n + j] -= L[i*n + k] * U[k*n + j];
			}
		}

        // cerr << "Printing the third A" << endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                // cerr << a[i * n + j] << " ";
//            }
//            // cerr << endl;
//        }
	}
}

int main (int argc, char* argv[]) {

	int n = atoi(argv[1]);
	int t = atoi(argv[2]);

	int i, j;

    double* a = new double[n * n];
    double* L = new double[n * n];
    double* U = new double[n * n];
    int* pi = new int[n];


    srand48(time(NULL));
    random_device rd;
//    mt19937 gen(rd());
    mt19937 gen(435);

    uniform_real_distribution<double> dis(0.0, 1.0);


    // generate random matrix
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// NOTE : replace with drand_r in the parallel versions
			// does this not require seeding?
            a[i*n + j] = dis(gen);
            a[i*n + j] = drand48();
			cout << a[i*n + j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";

    auto start = high_resolution_clock::now();
    decompose(a, pi, L, U, n);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<nanoseconds>(stop - start);
    cerr << duration.count() << endl;

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


	return 0;
}