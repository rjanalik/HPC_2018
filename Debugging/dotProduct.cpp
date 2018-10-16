#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#define NMAX 100
#define NMAX_ERR 1000000

using namespace std;

int main()
{
   double *a,*b;

   // Allocate memory for the vectors as 1-D arrays
   a = new double[NMAX];
   b = new double[NMAX];

   // Initialize the vectors with some values
   for(int i=0; i<NMAX; i++)
   {
      a[i] = i;
      b[i] = i/10.0;
   }

   // Compute dot product
   long double dotProduct = 0;
   for(int i=0; i< NMAX_ERR; i ++)
   {
      dotProduct += a[i] * b[i];
   }

   cout << "Dot product = " << dotProduct << endl;

   // De-allocate memory
   delete [] a;
   delete [] b;

   return 0;
}
