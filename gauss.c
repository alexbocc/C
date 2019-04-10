#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "lib.h"

int main(int argc, char *argv[])
{
    int i,n,m,j;
    int sorlen;
    mymat id,res,res2;
    unsigned char *words[100];
    char buff[1000];
    FILE  *fp;
    mymat a,c;
    mymat b,eq;
    if(argc!=2)
    {
        printf("Hasznalat:./matmul matrixFile vektorFile\n");
        exit(1);
    }
    a=readmat(argv[1],&n,&m);  
    id=ident(n);
    eq =appmat(a,id,n,m,n);
    gaussj(eq,n,n+m);
    printf("\n");
    printmat(eq,n,n+m);
    res = submat(eq,0,n-1,n,m+m-1);
    printf("\n");
    printmat(res,n,n);
    res2 = MatMatMul(res,n,n,a,n);
    printf("\n");
    printmat(res2,n,n);
    res2 = MatMatMul(a,n,n,res,n);
    printf("\n");
    printmat(res2,n,n);
    return 0;
}


