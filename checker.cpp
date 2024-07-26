#include "iostream"
#include "math.h"
#include <fstream>

using namespace std;

int main (int argc, char* argv[]) {
	
    ifstream file("dump.txt");
	int n, i, j, k;
	file >> n;
	int pi[n];
	double a[n*n], L[n*n], U[n*n];

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            file >> a[i*n + j];
        }
    }

	for (i = 0; i < n; i++) {
		 file >> pi[i];
	}



	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			file >> L[i*n + j];
		}
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			file >> U[i*n + j];
		}
	}

	double res[n*n];

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			res[i*n + j] = a[pi[i]*n + j];
			for (k = 0; k < n; k++) {
				res[i*n + j] -= L[i*n + k]*U[k*n + j];
			}
		}
	}

	// compute the L2,1 norm
	double norm = 0;
	for (j = 0; j < n; j++) {
		double col = 0;
		for (i = 0; i < n; i++) {
			col += res[i*n + j] * res[i*n + j];
		}
		norm += sqrt(col);
	}

	cout << norm;

	return 0;
}