#include<iostream>
//#include<conio.h>
using namespace std;

double largest(double arr[], int n)
{
    int i;

    // Initialize maximum element
    double max = arr[0];

    // Traverse array elements
    // from second and compare
    // every element with current max
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];

    return max;
}

int main(void) {
   double a[100][100], b[100], m[100], n[100], err[100];
   int p = 0, q = 0, i = 0, j = 0;
   cout << "Enter size of 2D array : ";
   cin >> p;
   for (i = 0; i < p; i++) {
      for (j = 0; j < p; j++) {
         cout << "a[" << i << ", " << j << " ]=";
         cin >> a[i][j];
      }
   }
   cout << "\nEnter values to the right side of equation\n";
   for (i = 0; i < p; i++) {
      cout << "b[" << i << ", " << j << " ]=";
      cin >> b[i];
   }
   cout << "Enter initial values of x\n";
   for (i = 0; i < p; i++) {
      cout << "x:[" << i<<"]=";
      cin >> n[i];
      m[i] = n[i]+1;
      err[i] = m[i]-n[i];
   }

   //cout << "\nEnter the no. of iteration : ";
   //cin >> q;
   double tol = 1e-6;
   while (largest(err, p)> tol) {
     //cout<<""
      for (i = 0; i < p; i++) {
         m[i] = n[i];
         n[i] = (b[i] / a[i][i]);
         for (j = 0; j < p; j++) {
            if (j == i)
               continue;
            n[i] = n[i] - ((a[i][j] / a[i][i]) * m[j]);
            err[i] = m[i]-n[i];
         }
         cout<<"x"<<i + 1 << "="<<n[i]<<" ";
      }
      cout << "\n";
      q--;
   }
   return 0;
}
