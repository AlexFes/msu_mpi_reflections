#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "reflections_MPI.cpp"
#include "reflections_MPI.h"

int main (int argc, char *argv[]) {

	MPI_Init (&argc, &argv);

	int t, T;

    int n, m, result = 1, x = 0;
	double time, discrepancy;

	char *input = NULL;

	MPI_Comm_rank (MPI_COMM_WORLD, &t);
	
	MPI_Comm_size (MPI_COMM_WORLD, &T);

	if ((argc!=3 && argc!=4) || !(n = atoi(argv[1])) || !(m = atoi(argv[2]))) {

        printf("\nВведите количество уравнений и размер блока.\n");
		return 0;
	}

	if (argc==4) {
	
		input = argv[3];
		
		x = 1;
	}

	int l = n%m;
	int N = n/m;
	
	if (l!=0) N++;



    double *a = new double [(n*n)/T + n*m];		//для системы
    double *b = new double [n + n*m];

	double *d = new double [m*(m+1)*N];			//для матриц отражения

	double *c = new double [m*m];				//для блоков
	double *c1 = new double [m*m];
	double *c2 = new double [m*m];
	

    readMatrix (n, m, t, T, a, b, input, x);

    if (n <= 20) writeMatrix (n, m, t, T, a, b, result);



	MPI_Barrier (MPI_COMM_WORLD);
	
	time = MPI_Wtime ();

    result = Function (n, m, t, T, a, b, d, c, c1, c2);

	time = MPI_Wtime () - time;



    if (n <= 20) writeMatrix (n, m, t, T, a, b, result);

	readMatrix (n, m, t, T, a, d, input, x);

	discrepancy = Discrepancy (a, b, t, T, d, c, c1, c2, b + n, n, m);



    if (t == 0 && result == 1) {
	
		printf ("\nВремя работы: %.2f\n", time);

        printf ("\nНорма невязки: %e\n\n", discrepancy);
	}

	else if (t == 0 && result == 0) 

		printf ("\nАлгоритм неприменим.\n\n");



	delete[]a;
	delete[]b;
	delete[]d;
	delete[]c;
	delete[]c1;
	delete[]c2;

	MPI_Finalize();

	return 0;
}
