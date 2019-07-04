#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>

using namespace std;


int** randGenArray(int rows, int cols){
   int res[rows][cols];
   for(int i = 0; i < rows; i++){
       for(int j = 0; j < cols; j++){
           res[i][j] = static_cast<int>(rand() % 5000);
       }
   } 

   return res;
}

int main(){

    int h = 100;
    int j = 100;
    int k = 100;

    int* a = randGenArray(h,j);
    int* b = randGenArray(j,k);
    int* res = randGenArray(h,k);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 100, 100, 100, 1, a, 100, b, 100, 0, res, 100);

    
}