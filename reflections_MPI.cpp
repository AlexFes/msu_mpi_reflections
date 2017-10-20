#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "reflections_MPI.h"

#define eps 1e-16

int max (int a, int b) {

	if (a > b) return a;
	else return b;
}	

double getElement (int n, int m, int t, int T, double *a, int i_local, int j_local) {			//получаем элемент	

	int p_local = j_local/m;
	int q_local = i_local/m;
	int l = n%m;

	int p_global = t + T*p_local;

	int Aj = (m*(p_global + 1) > n ? l : m);
	
return a[p_local*m*n + q_local*m*Aj + (i_local-q_local*m)*Aj + (j_local-p_local*m)];		
}												//legit

void getBlock (double *a, double *c, int Num, int Num1, int Num2) {			//получаем блок

	for (int i = 0; i < Num1; ++i) 

		for (int j = 0; j < Num2; ++j) {	
			c[i*Num2 + j] = a[Num];
			++Num;
		}

return;
}																	//legit

void putBlock (double *a, double *c, int Num, int Num1, int Num2) {

	for (int i = 0; i < Num1; ++i) 

		for (int j = 0; j < Num2; ++j) {	
			a[Num] = c[i*Num2 + j];
			++Num;
		}

return;
}

void readMatrix (int n, int m, int t, int T, double *a, double *b, char *input, int x) {			//считываем систему

	int Ai, Aj, bi=0, bj=t*m;
    MPI_Status st;

	int N_local, N_global = n/m;
	int l = n%m;

	if (l!=0) ++N_global;
	
	if (N_global%T != 0) N_local = N_global/T + (t < N_global%T ? 1 : 0);
	else N_local = N_global/T;

	FILE *in = fopen (input, "r");

	if (x == 0) {															//задаем матрицу по формуле
	
		for (int q = 0; q < N_global; ++q) {
	
			Ai = (m*(q+1) > n ? l : m);

			for (int ai = 0; ai < Ai; ++ai) {

				for (int p = 0; p < N_local; ++p) {

					Aj = (m*((t + T*p) + 1) > n ? l : m);

					for (int aj = p*m; aj < p*m + Aj; ++aj) {

						a[p*m*n + q*m*Aj + ai*Aj + (aj-p*m)] = 
                    //    1.0/(bi + bj + 1);
                    	n - max(bi,bj);
                    //  	abs(bi-bj);
						bj++;
					}

					bj+=(T-1)*m;
				}
				
				bj = t*m;
						
				bi++;
			}
		}
		
		if (t == 0) 
		
			for (int i = 0; i < n; i++) {
			
				b[i] = 0;	
						
				for (int j = 0; j < n; j+=2) b[i] += 
										//	1.0/(i + j + 1);
											n - max(i,j);
										//	abs(i-j);
			}
	}

	else if (x == 1)														//считываем матрицу из файла

		for (int q = 0; q < N_global; ++q) {
	
			Ai = (m*(q + 1) > n ? l : m);

			for (int ai = 0; ai < Ai; ++ai) {
			
				bj = 0;
				
				if (t == 0) b[bi] = 0;

				for (int p = 0; p < N_global; ++p) {

					Aj = (m*(p + 1) > n ? l : m);

					if (t == 0 && p%T == t) 

						for (int aj = (p/T)*m; aj < (p/T)*m + Aj; ++aj) { 

							fscanf (in, "%lf", &a[(p/T)*m*n + q*m*Aj + ai*Aj + (aj-(p/T)*m)]);

							if (bj%2 == 0) b[bi] += a[(p/T)*m*n + q*m*Aj + ai*Aj + (aj-(p/T)*m)];
							
							bj++;
						}

					else if (t == 0 && p%T != t) {

						for (int aj = p*m; aj < p*m + Aj; ++aj) {

							fscanf (in, "%lf", &b[n + (aj - p*m)]);

							if (bj%2 == 0) b[bi] += b[n + (aj - p*m)];
							
							bj++;
						}

						MPI_Send (b + n, Aj, MPI_DOUBLE, p%T, 0, MPI_COMM_WORLD);
					}
						
					else if (t != 0 && p%T == t) 

						MPI_Recv (a + (p/T)*m*n + q*m*Aj + ai*Aj, Aj, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
				}
				
				bi++;
			}
		}

	if (input != NULL) fclose(in);

return;
}																//legit

void writeMatrix (int n, int m, int t, int T, double *a, double *b, int result) {		//выводим систему
 
    int p_local, Aj, k = 0;
    MPI_Status st;
	
	int N_global = n/m;
	int l = n%m;

	if (l!=0) ++N_global;

    if (t == 0) if (result) printf("\n\n");

    for (int i = 0; i < n; i++) {

        if (t == 0) if (result) printf("\n");

        for (int p_global = 0; p_global*m < n; p_global++) {

            Aj = (m*(p_global + 1) > n ? l : m);
            p_local = p_global/T;

            if (t == 0 && p_global%T == 0)	for (int j = p_local*m; j < p_local*m + Aj; ++j) {
                	
                if (result) printf("%10f ", getElement (n, m, t, T, a, i, j));
            }

            else if (t == 0 && p_global%T != 0) { 

                MPI_Recv (b + n, Aj, MPI_DOUBLE, p_global%T, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

                for (int j = n; j < n + Aj; ++j) if (result) printf("%10f ", b[j]);
            }
            
            else if (t != 0 && p_global%T == t) {
            
            	for (int j = p_local*m; j < p_local*m + Aj; ++j) b[j - p_local*m] = getElement (n, m, t, T, a, i, j);

                MPI_Send (b, Aj, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            } 
        }
           
        if (t == 0) {            	    
        
        	if (result) printf("     %10f", b[k]);
        	k++;
        }
    }

    if (t == 0) if (result) printf("\n\n");
    
return;
}																//legit

void multiplyBlock (double *a, double *c, double *b, double *c2, bool k, int Num0, int Num1, int Num2, int Num3) {

	int i, j, p, NUM1, NUM3, tmpNum = Num0;
	double tmp1, tmp2, tmp3, tmp4;

	if (!(Num1 % 2)) NUM1 = Num1;
	else if (Num1 % 2) NUM1 = Num1 - 1;

	if (!(Num3 % 2)) NUM3 = Num3;
	else if (Num3 % 2) NUM3 = Num3 - 1;

	for (i = 0; i < Num2; ++i) 			

		for (j = 0; j < Num3; ++j) {

			b[i*Num3 + j] = a[tmpNum];			

			++tmpNum;
		}

	for (i = 0; i < NUM1; i += 2) 

		for (j = 0; j < NUM3; j += 2) {

			tmp1 = tmp2 = tmp3 = tmp4 = 0;

			for (p = 0; p < Num2; ++p) {
				tmp1 += c[i*Num2 + p]*b[p*Num3 + j];
				tmp2 += c[i*Num2 + p]*b[p*Num3 + j+1];
				tmp3 += c[(i+1)*Num2 + p]*b[p*Num3 + j];
				tmp4 += c[(i+1)*Num2 + p]*b[p*Num3 + (j+1)];
			}

			if (k) {
				a[Num0 + i*Num3 + j] = tmp1;
				a[Num0 + i*Num3 + j+1] = tmp2;
				a[Num0 + (i+1)*Num3 + j] = tmp3;
				a[Num0 + (i+1)*Num3 + j+1] = tmp4;
			}

			else if (!k) { 
				c2[i*Num3 + j] = tmp1;
				c2[i*Num3 + j+1] = tmp2;
				c2[(i+1)*Num3 + j] = tmp3;
				c2[(i+1)*Num3 + j+1] = tmp4;
			}
		}

	for (i = NUM1; i < Num1; i++) 	//когда нечетное число строк

		for (j = 0; j < Num3; j++) {

			tmp1 = 0;

			for (p = 0; p < Num2; ++p) tmp1 += c[i*Num2 + p]*b[p*Num3 + j];

			if (k) a[Num0 + i*Num3 + j] = tmp1;

			else if (!k) c2[i*Num3 + j] = tmp1;
		}

	for (i = 0; i < Num1; i++) 	//нечетное число столбцов

		for (j = NUM3; j < Num3; j++) {

			tmp1 = 0;

			for (p = 0; p < Num2; ++p) tmp1 += c[i*Num2 + p]*b[p*Num3 + j];

			if (k) a[Num0 + i*Num3 + j] = tmp1;

			else if (!k) c2[i*Num3 + j] = tmp1;
		}

return;
}

int inverseBlock (int m, double *a, double *b, double *c, int Num) {	//обращаем блок

    int i, j, k, tmpNum = 0;
    double tmp;

    for (i = 0; i < m; ++i)

        for (j = 0; j < m; ++j) {

            if (i == j) b[i*m + j] = 1;
            else b[i*m + j] = 0;
            c[i*m + j] = a[Num];
            ++Num;
        }

    for (i = 0; i < m; ++i) {

        tmp = 0;

        for (k = i; k < m; ++k)

            if (fabs(c[k*m + i]) > tmp) { tmp = c[k*m + i]; tmpNum = k; }

        if (fabs(tmp)<1e-14) return 0;

        for (k = 0; k < m; ++k) {

            tmp = c[i*m + k];
            c[i*m + k] = c[tmpNum*m + k];
            c[tmpNum*m + k] = tmp;

            tmp = b[i*m + k];
            b[i*m + k] = b[tmpNum*m + k];
            b[tmpNum*m + k] = tmp;
        }

        for (k = i + 1; k < m; ++k) {

            tmp = c[k*m + i]/c[i*m + i];

            for (j = i+1; j < m; ++j) c[k*m + j] -= tmp*c[i*m + j];

            for (j = 0; j < m; ++j) b[k*m + j] -= tmp*b[i*m + j];
        }
    }

    for (i = m - 1; i >= 0; --i) {

        for (j = 0; j < m; ++j) b[i*m + j] /= c[i*m + i];

        for (k = i - 1; k >= 0; --k)

            for (j = 0; j < m; ++j) b[k*m + j] -= b[i*m + j]*c[k*m + i];
    }

return 1;
}

void addBlock (double *a, double *c, int Num, int Num1, int Num2) {		//складываем блоки

	for (int i = 0; i < Num1; i++) 

		for (int j = 0; j < Num2; j++) {

			a[Num] -= c[i*Num2 + j];
			Num++;
		}

return;
}

double normBlock (double *a, int Num, int Num1, int Num2) {	//норма

	int i, j;
	double tmp, max = 0;

	for (i = 0; i < Num1; ++i) {

		tmp = 0;

		for (j = 0; j < Num2; ++j) {
			tmp += fabs (a[Num]);
			Num++;
		}		

		if (tmp > max) max = tmp;

	}

return max;
}

double Discrepancy (double *a, double *b, int t, int T, double *b1, double *c, double *c1, double *c2, double *res, int n, int m) {

    int root = 0, count = n, p_global, N_local, N = n/m;

	int N_global = n/m;
	int l = n%m;

	if (l!=0) ++N;
	
	if (N%T != 0) N_local = N/T + (t < N%T ? 1 : 0);
	else N_local = N/T;
	
	if (l!=0 && (N-1)%T == t) N_local--;

	for (int i = 0; i < n; i++) res[i] = 0;

	MPI_Bcast (b, count, MPI_DOUBLE, root, MPI_COMM_WORLD);

	for (int q_global = 0; q_global < N_global; q_global++) {

		for (int p_local = 0; p_local < N_local; p_local++) {

			p_global = t + T*p_local;

			getBlock (a, c, q_global*(m*m) + p_local*(n*m), m, m);

			multiplyBlock (b, c, c1, c2, 0, p_global*m, m, m, 1);

			for (int i = 0; i < m; i++) res[q_global*m + i] += c2[i];
		}

		if (l!=0 && N_global%T == t) {

			getBlock (a, c, q_global*(m*l) + N_local*(n*m), m, l);

			multiplyBlock (b, c, c1, c2, 0, N_global*m, m, l, 1);

			for (int i = 0; i < m; i++) res[q_global*m + i] += c2[i];
		}
	}

	if (l!=0) {

		for (int p_local = 0; p_local < N_local; p_local++) {

			p_global = t + T*p_local;

			getBlock (a, c, N_global*(m*m) + p_local*(n*m), l, m);

			multiplyBlock (b, c, c1, c2, 0, p_global*m, l, m, 1);

			for (int i = 0; i < l; i++) res[N_global*m + i] += c2[i];
		}

		if (N_global%T == t) {

			getBlock (a, c, N_global*(m*l) + N_local*(n*m), l, l);

			multiplyBlock (b, c, c1, c2, 0, N_global*m, l, l, 1);

			for (int i = 0; i < l; i++) res[N_global*m + i] += c2[i];
		}
	}

	if (t == 0) MPI_Reduce (res, b, count, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);	
		
	else MPI_Reduce (res, 0, count, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

	if (t == 0) for (int i = 0; i < n; i++) b[i] -= b1[i];	

return (t == 0 ? normBlock (b, 0, n, 1) : 0);
}

void reflections (double *a, double *d, int m, int Num, int Num_d, int size) { 				//метод отражений для диагонального блока

	double s = 0, norm_a = 0, norm_x = 0;

	for (int i = 0; i < size-1; i++) {		//шаг алгоритма

		for (int p = i+1; p < size; p++) s += a[Num + p*size + i]*a[Num + p*size + i];
			
		norm_a = sqrt(a[Num + i*size + i]*a[Num + i*size + i] + s);
			
		norm_x = sqrt((a[Num + i*size + i]-norm_a)*(a[Num + i*size + i]-norm_a) + s);
			
        if (fabs(norm_x) > Eps) {
						
			d[Num_d + i*(m+1)] = (a[Num + i*size + i] - norm_a)/norm_x;				
						
			for (int p = i+1; p < size; p++) d[Num_d + i*(m+1) + (p-i)] = a[Num + p*size + i]/norm_x;
        }

        else for (int p = i; p < size; p++) d[Num_d + i*(m+1) + (p-i)] = 0;

		a[Num + i*size + i] = norm_a;

		for (int p = i+1; p < size; p++) a[Num + p*size + i] = 0;

		for (int q = i+1; q < size; q++) {	//умножаем на x
			
			s = 0;
						
			for (int p = i; p < size; p++) s += a[Num + p*size + q]*d[Num_d + i*(m+1) + (p-i)];
			
            if (fabs(s)>Eps) for (int p = i; p < size; p++) a[Num + p*size + q] = a[Num + p*size + q] - 2*s*d[Num_d + i*(m+1) + (p-i)];
		}
			
		norm_a = 0;	
		norm_x = 0;
		s = 0;
	}
		
	return;
} 														//legit
		
void multiply (double *a, double *d, double *c, int m, int Num_v, int Num, int Num_d, int size, int num_col, int chaka) {	//умножаем не диагональный блок на m векторов отражения
		
    double s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, norm_x, norm_a;
			
	int tmp_count, tmp;
	
	tmp_count = (chaka==1 ? Num_v-1 : Num_v);
    tmp = num_col%10;

    if (chaka==1)                                  		//первый блок в столбце
			
        for (int count = 0; count < tmp_count; count++) {

            for (int j = 9; j < num_col-tmp; j+=10) {

                s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = s9 = s10 = 0;

                for (int i = count+1; i < size+1; i++) {

                    s1 += a[Num + (i-1)*num_col + j]*d[Num_d + count*(m+1) + (i-1-count)];
                    s2 += a[Num + (i-1)*num_col + (j-1)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s3 += a[Num + (i-1)*num_col + (j-2)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s4 += a[Num + (i-1)*num_col + (j-3)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s5 += a[Num + (i-1)*num_col + (j-4)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s6 += a[Num + (i-1)*num_col + (j-5)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s7 += a[Num + (i-1)*num_col + (j-6)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s8 += a[Num + (i-1)*num_col + (j-7)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s9 += a[Num + (i-1)*num_col + (j-8)]*d[Num_d + count*(m+1) + (i-1-count)];
                    s10 += a[Num + (i-1)*num_col + (j-9)]*d[Num_d + count*(m+1) + (i-1-count)];
                }


                for (int i = count+1; i < size+1; i++) {

                    a[Num + (i-1)*num_col + j] = a[Num + (i-1)*num_col + j] - 2*s1*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-1)] = a[Num + (i-1)*num_col + (j-1)] - 2*s2*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-2)] = a[Num + (i-1)*num_col + (j-2)] - 2*s3*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-3)] = a[Num + (i-1)*num_col + (j-3)] - 2*s4*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-4)] = a[Num + (i-1)*num_col + (j-4)] - 2*s5*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-5)] = a[Num + (i-1)*num_col + (j-5)] - 2*s6*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-6)] = a[Num + (i-1)*num_col + (j-6)] - 2*s7*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-7)] = a[Num + (i-1)*num_col + (j-7)] - 2*s8*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-8)] = a[Num + (i-1)*num_col + (j-8)] - 2*s9*d[Num_d + count*(m+1) + (i-1-count)];
                    a[Num + (i-1)*num_col + (j-9)] = a[Num + (i-1)*num_col + (j-9)] - 2*s10*d[Num_d + count*(m+1) + (i-1-count)];
                }
            }

            for (int j = num_col-tmp; j < num_col; j++) {

                s1 = 0;

                for (int i = count+1; i < size+1; i++) s1 += a[Num + (i-1)*num_col + j]*d[Num_d + count*(m+1) + (i-1-count)];

				if (fabs(s1)>Eps)

                	for (int i = count+1; i < size+1; i++) a[Num + (i-1)*num_col + j] = a[Num + (i-1)*num_col + j] - 2*s1*d[Num_d + count*(m+1) + (i-1-count)];

            }
        }

	if (chaka==2)										//нижние блоки

		for (int count = 0; count < tmp_count; count++) {

            for (int j = 9; j < num_col-tmp; j+=10) {

                s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = s9 = s10 = 0;

                s10 += c[count*num_col + (j-9)]*d[Num_d + count*(m+1)];
                s9 += c[count*num_col + (j-8)]*d[Num_d + count*(m+1)];
                s8 += c[count*num_col + (j-7)]*d[Num_d + count*(m+1)];
                s7 += c[count*num_col + (j-6)]*d[Num_d + count*(m+1)];
                s6 += c[count*num_col + (j-5)]*d[Num_d + count*(m+1)];
                s5 += c[count*num_col + (j-4)]*d[Num_d + count*(m+1)];
                s4 += c[count*num_col + (j-3)]*d[Num_d + count*(m+1)];
                s3 += c[count*num_col + (j-2)]*d[Num_d + count*(m+1)];
                s2 += c[count*num_col + (j-1)]*d[Num_d + count*(m+1)];
                s1 += c[count*num_col + j]*d[Num_d + count*(m+1)];

				for (int i = 1; i < size+1; i++) {

                    s10 += a[Num + (i-1)*num_col + (j-9)]*d[Num_d + count*(m+1) + i];
                    s9 += a[Num + (i-1)*num_col + (j-8)]*d[Num_d + count*(m+1) + i];
                    s8 += a[Num + (i-1)*num_col + (j-7)]*d[Num_d + count*(m+1) + i];
                    s7 += a[Num + (i-1)*num_col + (j-6)]*d[Num_d + count*(m+1) + i];
                    s6 += a[Num + (i-1)*num_col + (j-5)]*d[Num_d + count*(m+1) + i];
                    s5 += a[Num + (i-1)*num_col + (j-4)]*d[Num_d + count*(m+1) + i];
                    s4 += a[Num + (i-1)*num_col + (j-3)]*d[Num_d + count*(m+1) + i];
                    s3 += a[Num + (i-1)*num_col + (j-2)]*d[Num_d + count*(m+1) + i];
                    s2 += a[Num + (i-1)*num_col + (j-1)]*d[Num_d + count*(m+1) + i];
                    s1 += a[Num + (i-1)*num_col + j]*d[Num_d + count*(m+1) + i];
                }
				
                c[count*num_col + (j-9)] = c[count*num_col + (j-9)] - 2*s10*d[Num_d + count*(m+1)];
                c[count*num_col + (j-8)] = c[count*num_col + (j-8)] - 2*s9*d[Num_d + count*(m+1)];
                c[count*num_col + (j-7)] = c[count*num_col + (j-7)] - 2*s8*d[Num_d + count*(m+1)];
                c[count*num_col + (j-6)] = c[count*num_col + (j-6)] - 2*s7*d[Num_d + count*(m+1)];
                c[count*num_col + (j-5)] = c[count*num_col + (j-5)] - 2*s6*d[Num_d + count*(m+1)];
                c[count*num_col + (j-4)] = c[count*num_col + (j-4)] - 2*s5*d[Num_d + count*(m+1)];
                c[count*num_col + (j-3)] = c[count*num_col + (j-3)] - 2*s4*d[Num_d + count*(m+1)];
                c[count*num_col + (j-2)] = c[count*num_col + (j-2)] - 2*s3*d[Num_d + count*(m+1)];
                c[count*num_col + (j-1)] = c[count*num_col + (j-1)] - 2*s2*d[Num_d + count*(m+1)];
                c[count*num_col + j] = c[count*num_col + j] - 2*s1*d[Num_d + count*(m+1)];

                for (int i = 1; i < size+1; i++) {
					
                    a[Num + (i-1)*num_col + (j-9)] = a[Num + (i-1)*num_col + (j-9)] - 2*s10*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-8)] = a[Num + (i-1)*num_col + (j-8)] - 2*s9*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-7)] = a[Num + (i-1)*num_col + (j-7)] - 2*s8*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-6)] = a[Num + (i-1)*num_col + (j-6)] - 2*s7*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-5)] = a[Num + (i-1)*num_col + (j-5)] - 2*s6*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-4)] = a[Num + (i-1)*num_col + (j-4)] - 2*s5*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-3)] = a[Num + (i-1)*num_col + (j-3)] - 2*s4*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-2)] = a[Num + (i-1)*num_col + (j-2)] - 2*s3*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + (j-1)] = a[Num + (i-1)*num_col + (j-1)] - 2*s2*d[Num_d + count*(m+1) + i];
                    a[Num + (i-1)*num_col + j] = a[Num + (i-1)*num_col + j] - 2*s1*d[Num_d + count*(m+1) + i];
                }

			}
			
            for (int j = num_col-tmp; j < num_col; j++) {
			
                s1 = 0;
	
                s1 += c[count*num_col + j]*d[Num_d + count*(m+1)];
				
                for (int i = 1; i < size+1; i++) s1 += a[Num + (i-1)*num_col + j]*d[Num_d + count*(m+1) + i];

				if (fabs(s1)>Eps){				

                	c[count*num_col + j] = c[count*num_col + j] - 2*s1*d[Num_d + count*(m+1)];
				
                	for (int i = 1; i < size+1; i++)  a[Num + (i-1)*num_col + j] =	a[Num + (i-1)*num_col + j] - 2*s1*d[Num_d + count*(m+1) + i]; 
                }
			}
		}

	if (chaka==0) 										//первый столбец

		for (int count = 0; count < tmp_count; count++) {
			
            s1 = norm_a = norm_x = 0;
			
            for (int p = 0; p < size; p++) s1 += a[Num + p*m + count]*a[Num + p*m + count];
			
            norm_a = sqrt(c[count*m + count]*c[count*m + count] + s1);
			
            norm_x = sqrt((c[count*m + count]-norm_a)*(c[count*m + count]-norm_a) + s1);
				
            if (fabs(norm_x)>eps) {
            
            	d[Num_d + count*(m+1)] = (c[count*m + count] - norm_a)/norm_x;
            
            	for (int p = 1; p < size+1; p++) d[Num_d + count*(m+1) + p] = a[Num + (p-1)*m + count]/norm_x;
            }

            else for (int p = 0; p < size+1; p++) d[Num_d + count*(m+1) + p] = 0;
            
            if (fabs(s1)>eps) {
            
            	c[count*num_col + count] = norm_a;
            
            	for (int i = 1; i < size+1; i++) a[Num + (i-1)*num_col + count] = 0;
            }
            
            for (int j = count+1; j < num_col; j++) {													//умножаем матрицу 
	
                s1 = 0;

                s1 += c[count*num_col + j]*d[Num_d + count*(m+1)];
				
                for (int i = 1; i < size+1; i++) s1 += a[Num + (i-1)*num_col + j]*d[Num_d + count*(m+1) + i];			//вычисляем (a,x)

                if (fabs(s1)>eps) {																						//умножаем вектор
				
                    c[count*num_col + j] = c[count*num_col + j] - 2*s1*d[Num_d + count*(m+1)];
				
                    for (int i = 1; i < size+1; i++)  a[Num + (i-1)*num_col + j] =	a[Num + (i-1)*num_col + j] - 2*s1*d[Num_d + count*(m+1) + i];
                }
			}
		}

	return;
}			//legit
	
int Function (int n, int m, int t, int T, double *a, double *b, double *d, double *c, double *c1, double *c2) {

	int root, count, p_local, N_local, N = n/m, result = 1;
	MPI_Status st;
	
	int N_global = n/m;
	int l = n%m;

	if (l!=0) ++N;
	
	if (N%T != 0) N_local = N/T + (t < N%T ? 1 : 0);
	else N_local = N/T;
	
	if (l!=0 && (N-1)%T == t) N_local--;

	for (int q_global = 0; q_global < N_global; q_global++) {							//шаг алгоритма (приведение к верхнетреугольному виду)

		if (q_global%T == t) { 					//столбец наш 

			p_local = q_global/T;

			reflections (a, d, m, q_global*(m*m) + p_local*(n*m), q_global*m*(m+1), m); 		//метод отражений для диагонального блока
		
			for (int row = q_global+1; row < N_global; row++) 							//находим матрицы отражения (преобразуем первый столбец)

	    		multiply (a, d, a + q_global*(m*m) + p_local*(n*m), m, m, row*(m*m) + p_local*(n*m), row*m*(m+1), m, m, 0);

        	if (l!=0) multiply (a, d, a + q_global*(m*m) + p_local*(n*m), m, m, N_global*(m*m) + p_local*(n*m), N_global*m*(m+1), l, m, 0);
        	
        	count = (N-q_global)*m*(m+1);
        	
        	root = t;
        	
        	MPI_Bcast (d + q_global*m*(m+1), count, MPI_DOUBLE, root, MPI_COMM_WORLD); 
        	
        	p_local++;       	
		}

		else {									//столбец не наш

			p_local = ((q_global/T)*T + t < q_global ? q_global/T + 1 : q_global/T);

        	count = (N-q_global)*m*(m+1);

			root = q_global%T;

			MPI_Bcast (d + q_global*m*(m+1), count, MPI_DOUBLE, root, MPI_COMM_WORLD);
		}
					
		for (int col = p_local; col < N_local; col++) {							//преобразуем столбцы
		
	    	multiply (a, d, c, m, m, q_global*(m*m) + col*(n*m), q_global*m*(m+1), m, m, 1);

			for (int row = q_global+1; row < N_global; row++)
		
		    	multiply (a, d, a + q_global*(m*m) + col*(n*m), m, m, row*(m*m) + col*(n*m), row*m*(m+1), m, m, 2);

			if (l!=0) multiply (a, d, a + q_global*(m*m) + col*(n*m), m, m, N_global*(m*m) + col*(n*m), N_global*m*(m+1), l, m, 2);
		}

		if (l!=0 && N_global%T == t) {

			multiply (a, d, c, m, m, q_global*(m*l) + N_local*(n*m), q_global*m*(m+1), m, l, 1);
		
			for (int row = q_global+1; row < N_global; row++) 							//последний столбец
		
	    		multiply (a, d, a + q_global*(m*l) + N_local*(n*m), m, m, row*(m*l) + N_local*(n*m), row*m*(m+1), m, l, 2);

        	multiply (a, d, a + q_global*(m*l) + N_local*(n*m), m, m, N_global*(m*l) + N_local*(n*m), N_global*m*(m+1), l, l, 2);        
		}

		if (t == 0) {

			multiply (b, d, c, m, m, q_global*m, q_global*m*(m+1), m, 1, 1);
		
			for (int row = q_global+1; row < N_global; row++) 								//преобразуем столбец B
		
	    		multiply (b, d, b + q_global*m, m, m, row*m, row*m*(m+1), m, 1, 2);

        	if (l!=0) multiply (b, d, b + q_global*m, m, m, N_global*m, N_global*m*(m+1), l, 1, 2);
        }
	}		

	if (l!=0 && N_global%T == t) {																//метод отражений для последнего остаточного блока
	
        reflections (a, d, m, N_global*(m*l) + N_local*(n*m), N_global*m*(m+1), l);
	
		if (t!=0 && l > 1) MPI_Send (d + N_global*m*(m+1), m*(m+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
		else if (t == 0 && l > 1) multiply (b, d, c, m, l, N_global*m, N_global*m*(m+1), l, 1, 1);
	}
		
    else if (l > 1 && N_global%T!=t && t == 0) {
    	
    	MPI_Recv (d + N_global*m*(m+1), m*(m+1), MPI_DOUBLE, N_global%T, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
    	
    	multiply (b, d, c, m, l, N_global*m, N_global*m*(m+1), l, 1, 1);
    }	
 
	if (l!=0 && t == 0) { 													//обратный ход метода Гаусса
 
 		if (N_global%T == 0) {
 
            if (!inverseBlock (l, a, c, c1, N_global*(m*l) + N_local*(n*m))) result = 0;
        
        	multiplyBlock (b, c, c1, c2, 1, N_global*m, l, l, 1);
        
        	for (int q = N_global-1; q >= 0; q--) { 
        
        		getBlock (a, c, q*(m*l) + N_local*(n*m), m, l); 
        
        		multiplyBlock (b, c, c1, c2, 0, N_global*m, m, l, 1);
        	
        		addBlock (b, c2, q*m, m, 1);
			}
		}

		else if (N_global%T!=0) {
		
			MPI_Recv (b + n, n*l, MPI_DOUBLE, N_global%T, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		
            if (!inverseBlock (l, b + n, c, c1, N_global*(m*l))) result = 0;

			multiplyBlock (b, c, c1, c2, 1, N_global*m, l, l, 1);		
		
        	for (int q = N_global-1; q >= 0; q--) { 
        
        		getBlock (b + n, c, q*(m*l), m, l); 
        
        		multiplyBlock (b, c, c1, c2, 0, N_global*m, m, l, 1);
        	
        		addBlock (b, c2, q*m, m, 1);
			}
		}
	}	
		
	else if (l!=0 && t == N_global%T) MPI_Send (a + N_local*(n*m), n*l, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);	
		
    for (int q = N_global-1; q >= 0; q--) {   

		if (t == 0 && q%T == 0) {
		
            if (!inverseBlock (m, a, c, c1, q*(m*m) + (q/T)*(n*m))) result = 0;

        	multiplyBlock (b, c, c1, c2, 1, q*m, m, m, 1);

        	for (int p = q-1; p >= 0; p--) {

				getBlock (a, c, p*(m*m) + (q/T)*(n*m), m, m);

            	multiplyBlock (b, c, c1, c2, 0, q*m, m, m, 1);

            	addBlock (b, c2, p*m, m, 1);
        	}
    	}

		else if (t == 0 && q%T!=0) {

			MPI_Recv (b + n, n*m, MPI_DOUBLE, q%T, MPI_ANY_TAG, MPI_COMM_WORLD, &st);	
				
            if (!inverseBlock (m, b + n, c, c1, q*(m*m))) result = 0;

        	multiplyBlock (b, c, c1, c2, 1, q*m, m, m, 1);

        	for (int p = q-1; p >= 0; p--) {

				getBlock (b + n, c, p*(m*m), m, m);

            	multiplyBlock (b, c, c1, c2, 0, q*m, m, m, 1);

            	addBlock (b, c2, p*m, m, 1);
        	}
    	}
			
		else if (t!=0 && q%T == t) MPI_Send (a + (q/T)*(n*m), n*m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}	

	return result;
}