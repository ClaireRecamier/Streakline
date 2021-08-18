#include <stdio.h>
#include <stdlib.h>

int main () {
  //float ranvals[10000];
  FILE *fpt = fopen("RadialOffsets/ran4.csv", "w+");

  srand(0.0);


  for (int i = 0; i < 5;i++){
    //ranvals[i] = rand();
    fprintf(fpt,"%d\n",rand());
  }
  fclose(fpt);
}
