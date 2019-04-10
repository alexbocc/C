#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "lib.h"



int strbreak(unsigned char *buff,unsigned char **words,int dmin,int dmax)
{
  int i,fl;
  unsigned char *p;
  i=0;
  fl=1;
  for(p=buff;*p!=0;p++)
  {
    if((*p>=dmin) && (*p<=dmax))
    {
      *p=0;
      fl=1;
    }
    else if(fl==1)
    {
      fl=0;
      words[i]=p;
      i++;
    }

  }
  return i;
}

myvect initvect(int len)
{
  myvect v;
  v=(double *)malloc(len*sizeof(double));
  return v;
}


mymat initmat(int N,int M)
{
  int i;
  double *p;
  double *buff;
  mymat sorok;
  sorok=(double **)malloc(N*sizeof(double *));
  buff=(double *)malloc(N*M*sizeof(double));
  for(i=0,p=buff;i<N;i++,p=p+M)
  {
    sorok[i]=p;
  }
  return sorok;
}

void freevect(myvect *v)
{
  free(*v);

}

void freemat(mymat *w)
{
  free((*w)[0]);
  free(*w);
}


myvect readvect(char *name,int *size)
{
  int n,k,l,j;
  unsigned char *words[MAXSIZE];
  myvect v;
  FILE  *fp;
  char *buffer;
  struct stat st;
  k=stat(name, &st);
  if(k<0)
  {
    return NULL;
  }
  buffer=(char *)malloc(st.st_size+1);
  fp=fopen(name,"r");
  if(fp==NULL)
  {
    return NULL;
  }
  if(fread(buffer,1,st.st_size+1,fp)!=0);
  n=strbreak(buffer,words,0,32);
  v=initvect(n);
  *size=n;
  if(v==NULL)
  {
    return NULL;
  }
  for(k=0;k<n;k++)
  {
    v[k]=atof(words[k]);
  }
  fclose(fp);
  free(buffer);
  return v;
}

mymat readmat(char *name,int *N,int *M)
{

  int n,k,l,j;
  unsigned char **words;
  FILE *fp;
  char *buffer;
  struct stat st;
  mymat mm;
  k=stat(name, &st);
  if(k<0)
  {
    return NULL;
  }
  buffer=(char *)malloc(st.st_size+1);
  words = (unsigned char **)malloc(st.st_size*sizeof(char *)/2+1);
  fp=fopen(name,"r");
  if(fp==NULL)
  {
    free(buffer);
    free(words);
    return NULL;

  }
  l=0;
  while(fgets(buffer,st.st_size+1,fp)!=NULL)
  {
    n=strbreak(buffer,words,0,32);
    if(n==0)
    {
      continue;
    }
    if(l==0)
    {
      *M=n;
    }
    else if(n!=*M)
    {
      fclose(fp);
      free(buffer);
      free(words);
      return NULL;
    }
    l++;
  }
  *N=l;
  mm=initmat(*N,*M);
  fseek(fp,0,SEEK_SET);
  l=0;
  while(fgets(buffer,st.st_size+1,fp)!=NULL)
  {
    n=strbreak((char *)buffer,words,0,32);
    if(n==0)
    {
      continue;
    }
    for(j=0;j<n;j++)
    {
      mm[l][j]=atof(words[j]);
    }
    l++;
  }

  fclose(fp);
  free(buffer);
  free(words);
  return mm;

}

void matrixmul(mymat a,myvect b,myvect c,int n,int m)
{
  int i,j;
  for(i=0;i<n;i++)
  {
    c[i]=0;
    for(j=0;j<m;j++)
    {
      c[i]=c[i]+a[i][j]*b[j];
    }
  }
}

mymat MatMatMul(mymat m1, int m1rnum, int m1cnum, mymat m2, int m2cnum)
{
  mymat res;
  int i,j,k;

  res = initmat(m1rnum,m2cnum);
  for (i = 0 ; i < m1rnum ; i++)
  {
    for (j = 0 ; j < m2cnum ; j++)
    {
      res[i][j] = 0;
      for (k = 0 ; k < m1cnum ; k++)
      {
        res[i][j] += m1[i][k]*m2[k][j];
      }
    }
  }
  return res;
}

mymat appvect(mymat matrix, myvect v, int n, int m)
{
  int i,j,k;
  mymat y;
  y = initmat(n,m+1);
  for (i = 0 ; i < n ; i++)
  {
    for (j = 0  ; j < m ; j++)
    {
      y[i][j] = matrix[i][j];
    }
  }
  for (i = 0 ; i < n ; i++)
  {
    y[i][n] = v[i];
  }
  return y;
}

mymat appmat(mymat m1, mymat m2, int n, int m, int k)
{
  int i,j;
  mymat m3;
  m3=initmat(n,m+k);
  for (i = 0 ; i < n ; i++)
  {
    for (j = 0  ; j < m ; j++)
    {
      m3[i][j] = m1[i][j];
    }
  }
  for (i = 0 ; i < n ; i++)
  {
    for(j = 0  ; j < k ; j++)
    {
      m3[i][m+j] = m2[i][j];
    }
  }
  return m3;
}

mymat ident(int n)
{
  int i,j;
  mymat a;
  a=initmat(n,n);
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      if(i==j)
      {
        a[i][j]=1;
      }
      else
      {
        a[i][j]=0;
      }
    }

  }
  return a;
}

void chrow(mymat m1, int n, int k, int colnum)
{
  int i;
  double t;
  for(i=0; i < colnum; i++)
  {
    t=m1[n][i];
    m1[n][i]=m1[k][i];
    m1[k][i]=t;
  }
}

int gaussj(mymat m, int n, int k)
{
  int i,j,l,iv;
  double d;
  for(i=0;i<n-1;i++)
  {
    for(j=n-1;j>=i+1;j--)
    {
      if ((m[j-1][i] == 0) && (fabs(m[j][i]/m[j-1][i]) > 1e10))
      {
        chrow(m,j-1,j,k);
      }
      else
      {
        d=m[j][i]/m[j-1][i];
        for(l=i;l<k;l++)
        {
          m[j][l] -= d*m[j-1][l];
        }
      }

    }
//    printmat(m,n,k);
//    printf("\n");
  }

  for(i=n-1;i>0;i--)
  {
    for(j=0;j<=i-1;j++)
    {
      if ((m[j+1][i] == 0) && (fabs(m[j][i]/m[j+1][i]) > 1e10))
      {
        chrow(m,j+1,j,k);
      }
      else
      {
        d=m[j][i]/m[j+1][i];
        for(l=0;l<k;l++)
        {
          m[j][l] -= d*m[j+1][l];
        }

      }
    }
//    printmat(m,n,k);
//    printf("\n");
 }
//  printmat(m,n,k);
//  printf("\n");
  for (i = 0 ; i < n ; i++){
    if ((d=m[i][i]) != 0)
    {
      for (j = 0 ; j < k ; j++)
      {
        m[i][j] /= d;
      }
    }
    else
    {
      return -1;
    }
  }


  return 0;
}

mymat submat(mymat mt, int r1, int r2, int c1, int c2)
{
  int i,j,N,M;
  mymat mr;
  mr = initmat((N = r2-r1+1), (M = c2-c1+1));
  for (i = r1 ; i <= r2 ; i++)
  {
    for (j = c1 ; j <= c2 ; j++)
    {
      mr[i-r1][j-c1] = mt[i][j];
    }
  }
  return mr;
}

myvect column(mymat mt, int n, int m)
{
  int i,j,N,M;
  myvect u;
  u = initvect(m);
  for (i=0; i<m; i++)
  {
    u[i]=mt[i][n];
  }
  return u;
}


void printmat(mymat mt, int n, int m)
{
  int i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<m;j++)
    {
      printf("%13.4e ",mt[i][j]);
    }
    printf("\n");
  }

}


double poli(int iv, int n, myvect x)
{
  int i,r,k;
  double l;

  if (iv == 0)
  {
    r = 0;
    k = 0;
  }
  else
  {
    r = (iv-1)/n+1;
    k = (iv-1) % n;
  }
  l=1;
  for(i=0;i<r;i++)
  {
    l=l*x[k];
  }

  return l;
}

mymat approxmat(mymat x, int n, int m, int rnd)
{
  int i,j,k,iv,l,matsize;
  mymat ret;
  myvect v,u;
  v=initvect(m-2);
  u=initvect(m-2);
  matsize=1+rnd*(m-2);
  ret=initmat(matsize,matsize);


  for(i=0;i<matsize;i++)
  {
    for(j=0;j<matsize;j++)
    {
      ret[i][j] = 0;
      for (l = 0 ; l < n ; l++)
      {
        for(k=0;k<m-2;k++)
        {
          u[k]=x[l][k];
        }
        ret[i][j]+=poli(i,m-2,u)*poli(j,m-2,u);
      }
    }
  }


  return ret;
}

myvect approxvect(mymat x, int n, int m, int rnd)
{
  int i,j,k,l,vectsize;
  myvect y,u,v;
  vectsize=1+rnd*(m-2);
  y=initvect(vectsize);
  u=initvect(n);
  v=initvect(n);
  for(l=0;l<n;l++)
  {
    v[l]=x[l][m-2]+x[l][m-1];
  }
  for(i=0;i<vectsize;i++)
  {
    y[i]=0;
    for(j=0;j<n;j++)
    {
      for(k=0;k<m-2;k++)
      {
        u[k]=x[j][k];
      }
      y[i]=y[i]+poli(i,m-2,u)*v[j];
    }
  }
  return y;
}

double Newton(double x1, double x2)
{
  return ((GCONSTANT)*(EARTHMASS))/(x1*x1+x2*x2);
}
double NewtonX1(double x1, double x2)
{
  return -Newton(x1,x2)*x1/sqrt(x1*x1+x2*x2);
}
double NewtonX2(double x1, double x2)
{
  return -Newton(x1,x2)*x2/sqrt(x1*x1+x2*x2);;
}

void eulernextV(double x1, double x2, double v1, double v2, double *u1, double *u2, double dt)
{
  *u1=v1+NewtonX1(x1,x2)*dt;
  *u2=v2+NewtonX2(x1,x2)*dt;
}

void eulernextX(double x1, double x2, double v1, double v2, double *y1, double *y2, double dt)
{
  *y1=x1+v1*dt;
  *y2=x2+v2*dt;
}

double *RungeKutta(double *y1, double *y0, double t0, int dim, double t, double h, int (*jo)(double, double *, int ))
{
  int i;
  for (i = 0; i < dim; i++)
  {
    y1[i] = y0[i] + h * (jo(t0, y0, dim) / 6) + h * (jo(t0 + 0.5*h, y0 + 0.5*h*jo(t0, y0,dim),dim) / 3 + h * (jo(t0 + 0.5*h, y0 + 0.5*h*jo(t0 + 0.5*h, y0 + 0.5*h*jo(t0, y0,dim),dim),dim)) / 3 + (h / 6)*jo(t0 + h, y0 + h * jo(t0 + 0.5*h, y0 + 0.5*h*jo(t0 + 0.5*h, y0 + 0.5*h*jo(t0, y0))));
  }
  return NULL;
  

int jo(double t, double *y, int n)
{
  double r = 0;
  return r; 
}


