#include <stdio.h>
#include <stdlib.h>
#include "lib.h"

int main()
{
  int i;
  double x10,x20,v10,v20;
  double x1,x2,v1,v2,t,dt;
  FILE *f1;
  x10=405500000;
  x20=0;
  v10=0;
  v20=964;
  t=0;
  dt=1;
  f1=fopen("moonOrbit.txt","w+");
//  for(i=0;i<40320;i++)
  for(i=0;i<3600*24*29;i++)
  {
    eulernextV(x10,x20,v10,v20,&v1,&v2,dt);
    eulernextX(x10,x20,(v1+v10)/2,(v2+v20)/2,&x1,&x2,dt);
    fprintf(f1,"%12.1f %12.1f %12.5f\n",x1,x2,t/86400);
    x10=x1;
    x20=x2;
    v10=v1;
    v20=v2;
    t=t+dt;  
  }
  fclose(f1);

  return 0;
}
