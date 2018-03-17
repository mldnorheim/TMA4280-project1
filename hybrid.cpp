#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <numeric>
#include <functional>


using namespace std;

int main(int argc, char **argv)   // RIEMANN ZETA:

/* You can add an OpenMP directive to the same MPI program. The OpenMP parallelizeation is only active if you pass -fopenmp to
the compiler and will not affect the MPI execution, they are independent. */
{
    if (argc < 3){
        printf("Need more arguments: Remembered nthreads and n?\n");
        return 1;
    }

    int nprocs, myrank;
    const double pi = 3.141592653589793;
    int nthreads = atoi(argv[1]);
   // int n = atoi(argv[2]);

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


            double timestartzeta;
            int n;

            if (myrank == 0) {
                //cout << "Please enter a value for the number of approx. terms, n: " << endl;
                //cin >> n;
                n = atoi(argv[2]);
                timestartzeta = MPI_Wtime();
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


            double vR = 0.0;
            #pragma omp parallel for num_threads(nthreads) reduction(+:vR) //private(myid)
                for(int j=1; j<=n; j++){
                    vR += 1.0/((double)j*(double)j);  // sumR
                }

            double piR = sqrt(6*vR);
            double errorR;

            if (myrank==0) {
                double durationR = MPI_Wtime() - timestartzeta;
                double errorR = fabs(piR-pi);
                //cout << "MPI-Riemann approximates pi as: " << scientific << piR << endl;
                //cout << "The error is: " << scientific << errorR << endl;
                printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, errorR, durationR);
            }
   // MPI_Finalize();

    ////////////////////////////////////////////////////////////////////////////////////7


   /*
    int nprocs, myrank;
    const double pi = 3.141592653589793;
    int nthreads = atoi(argv[1]);

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   */

            double timestartmach;

            if (myrank == 0) {
                //cout << "Please enter a value for the number of approx. terms, n: " << endl;
                //cin >> n;
                n = atoi(argv[2]);
                timestartmach = MPI_Wtime();
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            double vM5 = 0.0;
            double vM239 = 0.0;
            #pragma omp parallel for num_threads(nthreads) reduction(+:vM5,vM239)
                for(int j=1; j<=n; j++){
                    // x=1/5
                    vM5 += pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);
                    // x=1/239
                    vM239 += pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);
                }

            double piM = 4.0*(4.0*vM5 - vM239);

            double errorM;
            if (myrank==0) {
                double durationM = MPI_Wtime() - timestartmach;
                double errorM = fabs(piM-pi);
                printf(" piM=%e \n errorM=%e \n durationM=%e \n", piM, errorM, durationM);
                //cout << "Machin approximates pi as: " << fixed << scientific << setprecision(15) << piM << endl;
                //cout << "The error is: " << scientific << errorM << endl;
            }

    MPI_Finalize();










    /////////////////////////////////////////////////////////////////////////////////////////

    return 0;

}
