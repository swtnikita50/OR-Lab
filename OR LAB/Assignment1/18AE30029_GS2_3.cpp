#include<iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void permute(string a, int l, int r)
{
    // Base case
    ofstream outfile;
   outfile.open("Permute.dat");
    if (l == r)

        outfile<<a<<endl;
    else
    {
        // Permutations made
        for (int i = l; i <= r; i++)
        {

            // Swapping done
            swap(a[l], a[i]);

            // Recursion called
            permute(a, l+1, r);

            //backtrack
            swap(a[l], a[i]);
        }
    }
}

//GuassSeidel Fucntion


double GS(double** a, double* b, double* k, int n){

  double x[n], x0[n];
  for (int i = 0; i < n; i++) {
    x[i] = k[i]+1.0; x0[i] = k[i];
  }

  double tol = 1e-6;

  while (abs(x0-x)> tol) {
    for (int l = 0; l<n; l++){
    x0[l] = k[l];}
     for (int i = 0; i < n; i++) {
        k[i] = (b[i] / a[i][i]);
        for (int j = 0; j < n; j++) {
           if (j == i)
              continue;
           k[i] = k[i] - ((a[i][j] / a[i][i]) * x0[j]);

           x[i] = k[i];
        }
        cout<<"x"<<i + 1 << "="<<k[i]<<" ";
     }
     cout << "\n";
  }
}


int main(void) {

    int n, m;
    cout << "Enter the number of variables:" << endl;
    cin >> n;
    cout << "Enter the number of equations:" << endl;
    cin >> m;

   double a[100][100], b[100], x0[100], x[100], k[100], A[100][100], B[100];
   int numSol =1;

   for (int i = 0; i<m; i++){
     numSol = numSol*(n-i);
   }
   //int p = 0, q = 0, i = 0, j = 0;

   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         cout << "A[" << i << ", " << j << " ]=";
         cin >> A[i][j];
      }
   }
   cout << "\nEnter values to the right side of equation\n";
   for (int i = 0; i < m; i++) {
      cout << "B[" << i << " ]=";
      cin >> B[i];
   }
   cout << "Enter initial values of x\n";
   for (int i = 0; i < n; i++) {
      cout << "x0:[" << i<<"]=";
      cin >> k[i];
   }
   for (int i = 0; i < n; i++) {
     x[i] = k[i]+1.0;
   }

   double tol = 1e-6;
   int idx[n-m];

   //numSol = nPm;
   string str, newstr;
   int var0idx[n-m], r;
   for (int i = 1; i<n; i++){
   str.push_back('i');
 }

 cout << str << endl;

 permute(str, 0, n-1);
 outfile.close();

 //take first n-m variables to be 0 from the permutations
 //i.e. remove the columns corresponding to the first n-m variables
infile.open("Permute.dat");

 for (int i = 0; i<numSol; i++){
   infile >> newstr;
   for (int j = 0; j = m; j++){
   newstr.pop_back();}
   for (int l = 0; l<n-m; l++){
   var0idx[l] = newstr.pop_back();}


   for (int p = 0; p<m; p++){
     b[p] = B[p];
     r=0;
     for (int q = 0; q<n; q++){

       for(int s = 0;s<n-m; s++){
         if (q == var0idx[s]) q=q+1;
       }
       a[p][r] = A[p][q];
       r++;
     }
   }

   gaussS(a,b, k, m);


 }


   return 0;
}
