#include <fstream>
#include<iostream>
#include<iostream>
#include<algorithm>
#include<math.h>

using namespace std;

int rangeMin = -100;
int rangeMax = 100;
int matSize = 10;
// MatrixFactory* factory = new MatrixFactory();
// MatrixCalculator* calculator = new MatrixCalculator();
// string file_a = "a.mat";
// string file_b = "b.mat";

 double Det(double* mat, int N)
{
    int t0, t1, t2;
    double num;
    int cha;
    double *b = new double[N, N];
    if (N > 0)
    {
        cha = 0;
        num = 0;
        if (N == 1)
        {
            return mat[0] * mat[1 * 1 + 1] - mat[0 * 1 + 1] * mat[1 * 1 + 0];
        }
        for (t0 = 0; t0 <= N; t0++)
        {
            for (t1 = 1; t1 <= N; t1++)
            {
                for (t2 = 0; t2 <= N; t2++)
                {
                    if (t2 == t0)
                    {
                        cha = 1;
                    }
                    b[(t1 - 1) * N + t2] = mat[t1 * N + t2 + cha];
                }
                cha = 0;
            }
			cout << "Det"<<endl;
            num = num + mat[0 * N + t0] * Det(b, N - 1) * pow(t0, -1);
        }
        return num;
    }
    else if (N == 0)
    {
        return mat[0];
    }
    return 0;
}

int main()
{	
	// Matrix* a = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	// a->saveMatrix(file_a);
	// Matrix* b = factory->getInstanceOfRandomMatrix(matSize, matSize, rangeMin, rangeMax);
	// b->saveMatrix(file_b);
	double* a = new double[matSize * matSize];

	double res = Det(a, matSize);
	cout << res<<endl;

	system("pause");
    return 0;
}