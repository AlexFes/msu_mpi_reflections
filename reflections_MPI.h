int max (int a, int b);
double getElement (int n, int m, double *a, int i, int j);
void getBlock (double *a, double *c, int Num, int Num1, int Num2);
void putBlock (double *a, double *c, int Num, int Num1, int Num2);
void readMatrix (int n, int m, double *a, double *b, char *input);
void writeMatrix (int n, int m, double *a, double *b);
void multiplyBlock (double *a, double *c, double *b, double *c2, bool k, int Num0, int Num1, int Num2, int Num3);
double normBlock (double *a, int Num, int Num1, int Num2);
void addBlock (double *a, double *c, int Num, int Num1, int Num2);
int inverseBlock (int m, double *a, double *b, double *c, int Num);
void multiplyBlock (double *a, double *c, double *b, double *c2, bool k, int Num0, int Num1, int Num2, int Num3);
double Discrepancy (double *a, double *b, double *b1, double *c, double *c1, double *c2, double *res, int n, int m);
void reflections (double *a, double *d, int n, int m, int Num, int size);
void multiply (double *a, double *d, int m, int Num_v, int Num, int size, int num_col, int chaka);
void Cache_Machine_On (int n, int m, double *a, double *b, int q, int p, int size, bool boo);
void Cache_Machine_Off (int n, int m, double *a, double *b, int q, int p, int size);
int Function (int n, int m, double *a, double *b, double *d, double *c, double *c1, double *c2);
#define Eps 1e-9
