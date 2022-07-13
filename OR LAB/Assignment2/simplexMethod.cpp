#include <iostream>
#include <math.h>

using namespace std;

int indexofSmallestElement(double array[], int size)
{
    int index = 0;

    for(int i = 1; i < size; i++)
    {
        if(array[i] < array[index])
            index = i;
    }

    return index;
}

int main(){
  int n,m;

  cout << "Enter no. of variables:" << endl;
  cin >> n;
  cout << "Enter no. of constraints:" << endl;
  cin >> m;

  double x(n), z(m), a(m+1, n+1), aLR(n+1);

  cout << "Enter simplex Tableau:" << endl;
  for (i=0;i++;i<m+1){
    for(j=0; j++; j<n+1){
      cin << a(j,i);
      if(i == m){
        aLR(j) = a(j,i);
      }
    }
  }


  for(j=0; j++; j<n+1){
      aLR(j) = a(j,m);
  }
  idxk = indexofSmallestElement(aLR, n+1);
  if(aLR(idxk)>0|| aLR(idxk)==0){
    optVal = a(n,m);
    //optSln = (0....0);
  }
  else{
    for(j=0; j++; j<m+1){
          aPC(j) = a(j,idxk);
      }

  }




}
