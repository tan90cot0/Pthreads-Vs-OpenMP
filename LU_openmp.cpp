#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iostream"
#include "omp.h"
#include "random"

#include "cassert"
#include "chrono"


#define DEBUG_MODE

using namespace std;
using namespace std::chrono;

int thread_count = 1;

void decompose(double* a, int* pi, double* L, double* U, int n) {
    int i, j, k;

#pragma omp parallel for shared(a, pi, L, U) private(i, j, k) schedule(static)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i > j)
                U[i * n + j] = 0; // 0's below the diagonal
            else if (i == j)
                L[i * n + j] = 1; // 1's on the diagonal
            else
                L[i * n + j] = 0; // 0's above the diagonal
        }
        pi[i] = i;    // identity permutation
    }

    for (k = 0; k < n; k++) {


        vector<pair<int, double>> maxes(thread_count, make_pair(0, 0.0));










#pragma omp parallel shared(a, maxes)
        {
            double max_local = 0.0; // Initialize max_local as 0 for each thread
            int max_index = -1;

            #pragma omp for
            for (int i = k; i < n; i++) {
                if (max_local < fabs(a[i * n + k])) {
                    max_local = fabs(a[i * n + k]);
                    max_index = i;
                }
            }

            // Store the maximum value and its index for each thread
            maxes[omp_get_thread_num()] = std::make_pair(max_index, max_local);
        }

        double max = 0;
        int max_idx = 0;
       for (int i = 0; i < thread_count; i++) {
            if (max < maxes[i].second) {
                max = maxes[i].second;
                max_idx = maxes[i].first;
            }
        }

#ifdef DEBUG_MODE
        double max2 = 0;
        int max_idx2 = 0;
        for (i = k; i < n; i++) {
            if (max2 < fabs(a[i*n + k])) {
                max2 = fabs(a[i*n + k]);
                max_idx2 = i;
            }
        }

        assert(max_idx == max_idx2);
        assert(max == max2);

        //    cerr << "Max: " << max << " Max2: " << max2 << endl;

#endif

//        //    cerr << std::abs(max) << endl;
        if (std::abs(max) < 1e-6) {// comparisons among float is bad, replace with |max| < eps ?
            //    cerr << "Singular Matrix\n";
        }

        swap(pi[k], pi[max_idx]);

        //    cerr << "Printing the first A"<< endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                //    cerr << a[i * n + j] << " ";
//            }
//            //    cerr << endl;
//        }

#pragma omp parallel for private(j)
        for (j = 0; j < n; j++) {    // hmm there must be an STL for this
            swap(a[k * n + j], a[max_idx * n + j]);
            if (j < k)
                swap(L[k * n + j], L[max_idx * n + j]);
        }

        //    cerr << "Printing the second A"<< endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                //    cerr << a[i * n + j] << " ";
//            }
//            //    cerr << endl;
//        }


        U[k * n + k] = a[k * n + k];



#pragma omp parallel for shared(L, U, a) private(i)
        for (i = k + 1; i < n; i++) {
            L[i * n + k] = a[i * n + k] / U[k * n + k];
            U[k * n + i] = a[k * n + i];
        }

#pragma omp parallel for shared(a) private(i, j)
        for (i = k + 1; i < n; i++) {
            for (j = k + 1; j < n; j++) {
                a[i * n + j] -= L[i * n + k] * U[k * n + j];
            }
        }

        //    cerr << "Printing the third A"<< endl;
//        for(int i = 0; i < n; i++){
//            for(int j = 0; j < n; j++){
//                //    cerr << a[i * n + j] << " ";
//            }
//            //    cerr << endl;
//        }

    }


}

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <matrix_size> <threads>" << endl;
        return 1;
    }

    random_device rd;
//    mt19937 gen(rd());
    mt19937 gen(435);

    srand48(time(NULL));

    uniform_real_distribution<double> dis(0.0, 1.0);

    int n = atoi(argv[1]);
    thread_count = atoi(argv[2]);

    omp_set_num_threads(thread_count);

    int i, j;

    double* a = new double[n * n];
    double* L = new double[n * n];
    double* U = new double[n * n];
    int* pi = new int[n];

    // generate random matrix
//#pragma omp parallel for shared(a) private(i, j)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
//            a[i * n + j] = dis(gen);
            a[i * n + j] = drand48();
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cout << a[i * n + j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";


    auto start = high_resolution_clock::now();
    decompose(a, pi, L, U, n);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<nanoseconds>(stop - start);
    cerr << duration.count() << endl;

    for (i = 0; i < n; i++) {
        cout << pi[i] << " ";
    }
    cout << "\n\n";

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cout << L[i * n + j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cout << U[i * n + j] << " ";
        }
        cout << "\n";
    }

    delete[] a;
    delete[] L;
    delete[] U;
    delete[] pi;

    return 0;
}
