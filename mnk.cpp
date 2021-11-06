#include <iostream>
#include <cmath>


double* gaus(double** matrix, double* b, int n);

int main() {
    const int N = 11;
    const int m = 2;
    double x[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    double y[11] = { 3, 87, 156, 210, 238, 252, 239, 211, 158, 90, -5 };

    //double *POWERX = new double[2 * m];
    std::cout << "POWERX:\n";
    for (int k = 0; k < 2 * m; k++)
    {
        POWERX[k] = 0;
        for (int i = 0; i < N; i++)
        {
            POWERX[k] += pow(x[i], k + 1);
        }
    }
    for (int k = 0; k < 2 * m; k++){
        std::cout << POWERX[k] << " ";
    }
    std::cout << "\n";

    double **SUMX = new double*[m + 1];
    for (int i = 0; i < m + 1; i++)
    {
        SUMX[i] = new double[m + 1];
    }

    for (int l = 0; l < m + 1; l++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            if (j + l)
            {
                SUMX[l][j] = POWERX[l + j - 1];
            }
            else
            {
                SUMX[l][j] = N;
            }
        }
    }

    double *PRAW = new double[m + 1];
    for (int l = 0; l < m + 1; l++)
    {
        PRAW[l] = 0;
        for (int i = 0; i < N; i++)
        {
            PRAW[l] += y[i] * pow(x[i], l);
        }
    }

    double * a = gaus(SUMX, PRAW, m + 1);
    double S2 = 0;
    for (int i = 0; i < N; i++)
    {
        double sum = y[i];
        for (int j = 0; j < m + 1; j++)
        {
            sum -= a[j] * pow(x[i], j);
        }
        S2 += pow(sum, 2);
    }
    S2 /= N - m - 1;
    double sigma = sqrt(S2);
    std::cout << "a: \n";
    for (int i = 0; i < m + 1; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << "\nsigma: " << sigma;


    return 0;
}

double* gaus(double** matrix, double* b, int n){
    int ior[n], l;
    for (int k = 0; k < n; ++k) {
        ior[k] = k;
    }
    for (int k = 0; k < n; ++k) {
        double akk = 0;
        int M, p;
        for (int i = k; i < n; ++i) {
            l = ior[i];
            if (fabs(matrix[l][k]) < akk) {
                continue;
            }
            M = l;
            p = i;
            akk = fabs(matrix[l][k]);
        }
        ior[p] = ior[k];
        ior[k] = M;
        double amain = matrix[M][k];
        if (amain == 0) {
            std::cout << "error";
            return nullptr;
        }
        for (int j = k; j < n; ++j) {
            matrix[M][j] = matrix[M][j] / amain;
        }
        b[M] = b[M] / amain;
        for (int i = k + 1; i < n; ++i) {
            l = ior[i];
            for (int j = k + 1; j < n; ++j) {
                matrix[l][j] = matrix[l][j] - matrix[l][k] * matrix[M][j];
            }
            b[l] = b[l] - matrix[l][k] * b[M];
        }
    }
    l = ior[n - 1];
    b[l] = b[l] / matrix[l][n - 1];
    if (matrix[l][n - 1] == 0) {
        std::cout << "error";
        return nullptr;
    }
    auto* x = new double[n];
    x[n - 1] = b[l];
    for (int k = n - 2; k > -1; --k) {
        l = ior[k];
        double s = 0;
        for (int j = k + 1; j < n; ++j) {
            s += matrix[l][j] * x[j];
        }
        x[k] = b[l] - s;
    }
    return x;
}
