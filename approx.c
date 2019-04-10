#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <conio.h>
#include "lib.h" 

int main(int argc, char *argv[])
{
  int i,n,m,j,rnd,t,vectsize;
  mymat P,a,A;
  myvect c,x,jo;
  double y;
  FILE *f1;
  if(argc!=3)
  {
    printf("Hasznalat:./approx matrixFile kozelitesRend\n");
    exit(1);
  }
  rnd=atoi(argv[2]);
  a=readmat(argv[1],&n,&m);
  x=initvect(m-2);
  vectsize=1+rnd*(m-2);
  c=initvect(vectsize);
  jo=initvect(vectsize);
  P = approxmat(a,n,m,rnd);
  jo=approxvect(a,n,m,rnd);
  A=appvect(P,jo,vectsize,vectsize);
  gaussj(A,vectsize,vectsize+1);
  for(i=0;i<vectsize;i++)
  {
    c[i]=A[i][vectsize];
  }
  f1=fopen("fittedData.txt","w+");
  for(i=0;i<n;i++)
  {
    for(j=0;j<m-2;j++)
    {
      x[j]=a[i][j];
    }
    y=0;
    for(j=0;j<vectsize;j++)
    {
      y=y+c[j]*poli(j,m-2,x);
    }
    fprintf(f1,"%f %7.2f\n",y,a[i][m-2]+a[i][m-1]);
  }
   fclose(f1);
   for(i=0;i<vectsize;i++)
   {
     printf("%f\n",c[i]);
   }
   freemat(&P);
   freemat(&A);
   freemat(&a);
   freevect(&c);
   freevect(&x);
   freevect(&jo);

 return 0;
  
}