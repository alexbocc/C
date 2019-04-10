#define MAXSIZE 1000

typedef double *myvect;

typedef double **mymat;

int strbreak(unsigned char *buff,unsigned char **words,int dmin,int dmax);

myvect initvect(int len);

mymat initmat(int N,int M);

myvect readvect(char *name,int *size);

mymat readmat(char *name,int *N,int *M);

void freevect(myvect *v);

void freemat(mymat *w);

void matrixmul(mymat a,myvect b,myvect c,int n,int m);

mymat appvect(mymat matrix, myvect v, int n, int m);

mymat appmat(mymat m1, mymat m2, int n, int m, int k);

mymat ident(int n);

void chrow(mymat m1, int n, int k, int colnum);

int gaussj(mymat m, int n, int k);

void printmat(mymat mt, int n, int m);

mymat submat(mymat mt, int r1, int r2, int c1, int c2);

mymat MatMatMul(mymat m1, int m1rnum, int m1cnum, mymat m2, int m2cnum);

myvect column(mymat mt, int n, int m);

double poli(int iv, int n, myvect x);

mymat approxmat(mymat x, int n, int m, int rnd);

myvect approxvect(mymat x, int n, int m, int rnd);
