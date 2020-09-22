#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <iostream>
using namespace std;
using namespace Eigen;
int main()
{

    VectorXd x(2), b(2);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.push_back(T(0, 0, 3));
    tripletList.push_back(T(1, 0, 2.5));
    tripletList.push_back(T(0, 1, -1));
    tripletList.push_back(T(1, 1, 10));
    b(0) = 1;
    b(1) = 2;
    SparseMatrix<double> A(2, 2);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
    lscg.compute(A);
    x = lscg.solve(b);
    std::cout << "#iterations:     " << lscg.iterations() << std::endl;
    std::cout << "estimated error: " << lscg.error() << std::endl;
    std::cout << x << endl;
    // update b, and solve again
    x = lscg.solve(b);

    MatrixXd matA(2, 2), matB(2, 1), sol;
    matA(0, 0) = 3;
    matA(1, 0) = 2.5;
    matA(0, 1) = -1;
    matA(1, 1) = 10;

    matB(0, 0) = 1;
    matB(1, 0) = 2;
    cout << matA << endl;

    sol = matA.colPivHouseholderQr().solve(matB);

    cout << sol << endl;

    cin.get();
    return 0;
}