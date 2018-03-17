#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <numeric>
#include <functional>
//#pragma omp parallel default(shared) private(beta,pi)
using namespace std;

int main(int argc, char **argv)  // compiling command: -fopenmp
{

    if (argc < 3){
        printf("Need more arguments: Remembered nthreads and n?\n");
        return 1;
    }

    int nthreads = atoi(argv[1]);
    int n = atoi(argv[2]);

    const double pi = 3.141592653589793;
    double time_start = omp_get_wtime();
   // int myid = omp_get_thread_num();
    double sumR = 0.0;

    #pragma omp parallel for num_threads(nthreads) reduction(+:sumR) //private(myid)
        for(int j=1; j<=n; j++){
            sumR += 1.0/((double)j*(double)j);
        }

    double piR = sqrt(6*sumR);
    double error = fabs(piR-pi);
    double duration = omp_get_wtime() - time_start;

    printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);

    return 0; //piR;

 /*
     int nthreads, myid;
    const double pi = 3.141592653589793;

    #pragma omp parallel shared(n) private(myid)
    {
        double time_start = omp_get_wtime();

        myid = omp_get_thread_num();
      //  #pragma omp single shared(nthreads) private(myid)
        int n;

        if (myid==0){
           // nthreads = omp_get_num_threads();
           // printf("the number of threads are %d\n", nthreads);
            cout << "This is master thread, please enter a value for n: " << endl;
            cin >> n;

        }

        #pragma omp for reduction(+:sumR)
        double vR = 0.0;
        for(int j=myrank+1; j<=n; j+=nprocs){
            vR += 1.0/((double)j*(double)j);
        }
        double sumR;
        MPI_Reduce(&vR, &sumR, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // double SR = accumulate(vR.begin(),vR.end(),0.0);

        double piR = sqrt(6*sumR);
        double error;

        if (myrank==0) {
            double duration = MPI_Wtime() - time_start;
            double error = fabs(piR-pi);
            //cout << "MPI-Riemann approximates pi as: " << scientific << piR << endl;
            //cout << "The error is: " << scientific << errorR << endl;
            printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);
        }


    }
     .
 */

}




