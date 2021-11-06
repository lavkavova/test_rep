#include <iostream>
#include <iomanip>
#include <cmath>

void Jcobian(double x1, double x2, double **J);
void Jcobian_M(double x1, double x2, double **J);
double* gaus(double** matrix, double* b, int n);
double* neutone(double x1, double x2);
double* neutone_M(double x1, double x2);
double f1(double x1, double x2);
double f2(double x1, double x2);
double df1dx1(double x1, double x2);
double df1dx2(double x1, double x2);
double df2dx1(double x1, double x2);
double df2dx2(double x1, double x2);


const double eps = 1e-9;
const double M_const = 0.001;
const int NIT = 110;


int main() {
    double x1 = 1.0;
    double x2 = 1.0;
    int n = 2;
    auto* otvet = neutone(x1, x2);
    auto* otvet_M = neutone_M(x1, x2);

    std::cout << "\n";
    for (int i = 0; i < n; i++)
        std::cout << otvet[i] << std::endl;

    std::cout << "\n";
    for (int i = 0; i < n; i++)
        std::cout << otvet_M[i] << std::endl;

    return 0;
}

double f1(double x1, double x2){
    return 2*pow(x1, 3) + pow(x2, 2)-1;
}

double f2(double x1, double x2){
    return x1*pow(x2, 3) - x2 -4;
}

double df1dx1(double x1, double x2){
    return 6*pow(x1, 2);
}

double df1dx2(double x1, double x2){
    return 2*x2;
}

double df2dx1(double x1, double x2){
    return pow(x2, 3);
}

double df2dx2(double x1, double x2){
    return 3*x1*pow(x2, 2)-1;
}


void Jcobian(double x1, double x2, double **J){
    J[0][0] = df1dx1(x1, x2);
    J[0][1] = df1dx2(x1, x2);
    J[1][0] = df2dx1(x1, x2);
    J[1][1] = df2dx2(x1, x2);
}

void Jcobian_M(double x1, double x2, double **J){
    J[0][0] = (f1(x1 + M_const * x1, x2) - f1(x1, x2)) / (M_const * x1);
    J[0][1] = (f1(x1, x2 + M_const * x2) - f1(x1, x2)) / (M_const * x2);
    J[1][0] = (f2(x1 + M_const * x1, x2) - f2(x1, x2)) / (M_const * x1);
    J[1][1] = (f2(x1, x2 + M_const * x2) - f2(x1, x2)) / (M_const * x2);
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

double* neutone(double x1, double x2){
    int k = 1;
    std::cout << std::left << std::setw(4) << "k" << std::setw(20) << "d1" << std::setw(20) << "d2"
              << std::setw(20) << "x1" << std::setw(20) << "x2" << std::endl;
    auto *F = new double[2];
    auto **J = new double*[2];
    for (int i = 0; i < 2; i++){
        J[i] = new double[2];
    }
    auto* dX = new double[2];
    double x1k, x2k;
    double d1, d2;
    double tmp;
    do{
        F[0] = -f1(x1, x2);
        F[1] = -f2(x1, x2);
        Jcobian(x1, x2, J);

        dX = gaus(J, F, 2);
        x1k = x1 + dX[0];
        x2k = x2 + dX[1];
        d1 = fabs(f1(x1, x2));
        tmp = fabs(f2(x1, x2));
        if (tmp > d1){
            d1 = tmp;
        }
        d2 = fabs(x1k - x1) / (x1k >= 1 ? x1k : 1);
        tmp = fabs(x2k - x2) / (x2k >= 1 ? x2k : 1);
        if (tmp > d2){
            d2 = tmp;
        }
        x1 = x1k;
        x2 = x2k;
        std::cout << std::left << std::setw(4) << k << std::setw(20) << d1 << std::setw(20) << d2 << std::setw(20)
                  << x1 << std::setw(20) << x2 << std::endl;
        if (k >= NIT){
            std::cout << "IER=2\n";
            system("pause");
            exit(2);
        }
        k++;
    } while (d1>eps && d2>eps);
    dX[0] = x1;
    dX[1] = x2;
    return dX ;
}

double* neutone_M(double x1, double x2){
    int k = 1;
    std::cout << std::left << std::setw(4) << "k" << std::setw(20) << "d1" << std::setw(20) << "d2"
              << std::setw(20) << "x1" << std::setw(20) << "x2" << std::endl;
    auto *F = new double[2];
    auto **J = new double*[2];
    for (int i = 0; i < 2; i++){
        J[i] = new double[2];
    }
    auto *dX = new double[2];
    double x1k, x2k;
    double d1, d2;
    double tmp;
    do{
        F[0] = -f1(x1, x2);
        F[1] = -f2(x1, x2);
        Jcobian_M(x1, x2, J);

        dX = gaus(J, F, 2);
        x1k = x1 + dX[0];
        x2k = x2 + dX[1];
        d1 = fabs(f1(x1, x2));
        tmp = fabs(f2(x1, x2));
        if (tmp > d1){
            d1 = tmp;
        }
        
        d2 = fabs(x1k - x1) / (x1k >= 1 ? x1k : 1);
        tmp = fabs(x2k - x2) / (x2k >= 1 ? x2k : 1);
        if (tmp > d2){
            d2 = tmp;
        }
        x1 = x1k;
        x2 = x2k;
        std::cout << std::left << std::setw(4) << k << std::setw(20) << d1 << std::setw(20) << d2 << std::setw(20)
                  << x1 << std::setw(20) << x2 << std::endl;
        if (k >= NIT){
            std::cout << "IER=2\n";
            system("pause");
            exit(2);
        }
        k++;
    } while (d1>eps && d2>eps);
    dX[0] = x1;
    dX[1] = x2;
    return dX ;
}
