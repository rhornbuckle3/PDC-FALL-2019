#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <omp.h>
int main(void){

float a[7];
//a=malloc(sizeof(float)*7);
int i;
//i=malloc(sizeof(int))
for(i=0;i<7;i++){
    a[i]=(rand()*clock())%50;
}
printf("%s","Input Array\n");
for(i=0;i<7;i++){
    printf("%d",(int)a[i]);
    printf("%s","\n");
}
for(i=1;i<7;i++){
    a[i]=a[i-1]+a[i];
}
printf("%s","Output Array\n");
for(i=0;i<7;i++){
    printf("%d",(int)a[i]);
    printf("%s","\n");
}
return 0;
}