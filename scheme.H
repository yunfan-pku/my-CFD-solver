#ifndef SCHEME_H
#define SCHEME_H
#include<list>
#include "baseMesh.H"
#include "boundary.H"
#include "equation.H"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
using namespace Eigen;
#include <cmath>
class solverManger;
const int D_X = 0;
const int D_Y = 1;
class scheme
{
public:
    solverManger* psm;
    meshBase& m;
    int is, ie;
    int js, je;
    scheme(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input) :m(m_input), is(is_input), ie(ie_input), js(js_input), je(je_input), psm(nullptr) {}
    virtual void execute() {}
    virtual void init() {}
    virtual double residual() = 0;
};
namespace elliptic {
    class pointGaussSeidel :public scheme
    {
    public:
        ellipticEuqation& eqn;
        pointGaussSeidel(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, ellipticEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input) {}
        double compute(int i, int j) {
            /*
            double DX = (m.x[i + 1] - m.x[i - 1]) / 2;
            double DX2I = 1.0 / (DX * DX);
            double DY = (m.y[j + 1] - m.y[j - 1]) / 2;
            double DY2I = 1.0 / (DY * DY);
            return 1.0 / 2.0 / (eqn.A * DX2I + DY2I) * (eqn.A * (m.value[i + 1][j] + m.value[i - 1][j]) * DX2I + (m.value[i][j + 1] + m.value[i][j - 1]) * DY2I);
            */
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return 1 / (eqn.A * (DXI1 * DXI3 + DXI2 * DXI3) + DYI1 * DYI3 + DYI2 * DYI3) * (eqn.A * (m.now()[i + 1][j] * DXI1 * DXI3 + m.now()[i - 1][j] * DXI2 * DXI3) + (m.now()[i][j + 1] * DYI1 * DYI3 + m.now()[i][j - 1] * DYI2 * DYI3));

        }
        void execute()
        {
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    m.now()[i][j] = compute(i, j);
                }
        }
        double residualLoc(int i, int j)
        {
            /*
            double DX = (m.x[i + 1] - m.x[i - 1]) / 2;
            double DX2I = 1.0 / (DX * DX);
            double DY = (m.y[j + 1] - m.y[j - 1]) / 2;
            double DY2I = 1.0 / (DY * DY);
            return fabs(eqn.A * (m.value[i + 1][j] - 2 * m.value[i][j] + m.value[i - 1][j]) * DX2I + (m.value[i][j + 1] - 2 * m.value[i][j] + m.value[i][j - 1]) * DY2I);
    */
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return fabs(eqn.A * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3 + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                        resi = temp;
                }
            return resi;
        }

    };

    class pointJacobi :public scheme
    {
    public:
        ellipticEuqation& eqn;
        pointJacobi(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, ellipticEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input) {}
        double compute(int i, int j) {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return 1 / (eqn.A * (DXI1 * DXI3 + DXI2 * DXI3) + DYI1 * DYI3 + DYI2 * DYI3) * (eqn.A * (m.old()[i + 1][j] * DXI1 * DXI3 + m.old()[i - 1][j] * DXI2 * DXI3) + (m.old()[i][j + 1] * DYI1 * DYI3 + m.old()[i][j - 1] * DYI2 * DYI3));
        }
        void execute()
        {
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    m.now()[i][j] = compute(i, j);
                }
        }
        double residualLoc(int i, int j)
        {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return fabs(eqn.A * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3 + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                        resi = temp;
                }
            return resi;
        }

    };



    template<int dirc>
    class lineJacobi :public scheme
    {
    public:
        ellipticEuqation& eqn;
        lineJacobi(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, ellipticEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input) {}
        double residualLoc(int i, int j)
        {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return fabs(eqn.A * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3 + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            //int si, sj;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                    {
                        resi = temp;
                        //si = i;
                       // sj = j;
                    }
                }
            //cout << "!!!!!!!!!!!!" << si << "," << sj << endl;
            return resi;
        }
        void execute();
    };

    template<>
    void lineJacobi<D_X>::execute()
    {
        int n_mtr = ie - is + 2;
        VectorXd x(n_mtr), matB(n_mtr);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        int i;
        for (int j = js;j < je;j++)
        {
            tripletList.clear();
            i = 0;
            if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(i, i, 1.0));
                tripletList.push_back(T(i, i + 1, -1.0));
                /*matA(i, i) = 1.0;
                matA(i, i + 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(i) = -(m.x()[i + 1] - m.x()[i]) * np->f(m.x()[i], m.y()[j]);

            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(i, i, 1.0));
                /*
                 matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(i) = dp->f(m.x()[i], m.y()[j]);
            }
            for (i++;i < ie;i++)
            {
                double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
                double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
                double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
                double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
                double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
                double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
                tripletList.push_back(T(i, i, eqn.A * (DXI1 + DXI2) * DXI3 + (DYI1 + DYI2) * DYI3));
                tripletList.push_back(T(i, i - 1, -eqn.A * DXI2 * DXI3));
                tripletList.push_back(T(i, i + 1, -eqn.A * DXI1 * DXI3));
                /*
                matA(i, i) = eqn.A * (DXI1 + DXI2) * DXI2 + (DYI1 + DYI2) * DYI3;
                matA(i, i - 1) = -eqn.A * DXI2 * DXI2;
                matA(i, i + 1) = -eqn.A * DXI1 * DXI2;
                */
                matB(i) = (m.old()[i][j + 1] * DYI1 + m.old()[i][j - 1] * DYI2) * DYI3;

            }
            if (m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(i, i, 1.0));
                tripletList.push_back(T(i, i - 1, -1.0));
                /*
                matA(i, i) = 1.0;
                matA(i, i - 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(i) = (m.x()[i] - m.x()[i - 1]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(i, i, 1.0));
                /*
                matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(i) = dp->f(m.x()[i], m.y()[j]);
            }
            SparseMatrix<double> matA(n_mtr, n_mtr);
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
            lscg.compute(matA);
            x = lscg.solve(matB);
            */

            BiCGSTAB<SparseMatrix<double> > solver;
            solver.compute(matA);
            x = solver.solve(matB);


            //MatrixXd sol = matA.colPivHouseholderQr().solve(matB);
            /*
            SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
            solver.analyzePattern(matA);
            solver.factorize(matA);
            x = solver.solve(matB);
            */

            for (i = 0;i <= ie;i++)
                m.now()[i][j] = x(i);



        }
    }

    template<>
    void lineJacobi<D_Y>::execute()
    {
        int n_mtr = je - js + 2;
        VectorXd x(n_mtr), matB(n_mtr);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        int j;
        for (int i = is;i < ie;i++)
        {
            tripletList.clear();
            j = 0;
            if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j + 1, -1.0));
                /*matA(i, i) = 1.0;
                matA(i, i + 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = -(m.y()[j + 1] - m.y()[j]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                 matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            for (j++;j < je;j++)
            {
                double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
                double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
                double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
                double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
                double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
                double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
                tripletList.push_back(T(j, j, eqn.A * (DXI1 + DXI2) * DXI3 + (DYI1 + DYI2) * DYI3));
                tripletList.push_back(T(j, j - 1, -DYI2 * DYI3));
                tripletList.push_back(T(j, j + 1, -DYI1 * DYI3));
                /*
                matA(i, i) = eqn.A * (DXI1 + DXI2) * DXI2 + (DYI1 + DYI2) * DYI3;
                matA(i, i - 1) = -eqn.A * DXI2 * DXI2;
                matA(i, i + 1) = -eqn.A * DXI1 * DXI2;
                */
                matB(j) = eqn.A * (m.old()[i + 1][j] * DXI1 + m.old()[i - 1][j] * DXI2) * DXI3;

            }
            if (m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j - 1, -1.0));
                /*
                matA(i, i) = 1.0;
                matA(i, i - 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = (m.y()[i] - m.y()[i - 1]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            SparseMatrix<double> matA(n_mtr, n_mtr);
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
            lscg.compute(matA);
            x = lscg.solve(matB);
            */

            BiCGSTAB<SparseMatrix<double> > solver;
            solver.compute(matA);
            x = solver.solve(matB);

            //MatrixXd sol = matA.colPivHouseholderQr().solve(matB);
            /*
            SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
            solver.analyzePattern(matA);
            solver.factorize(matA);
            x = solver.solve(matB);
            */

            for (j = 0;j <= je;j++)
                m.now()[i][j] = x(j);
        }
    }


    template<int dirc>
    class gaussSeidelLineRelaxation :public scheme
    {
    public:
        ellipticEuqation& eqn;
        gaussSeidelLineRelaxation(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, ellipticEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input) {}
        double residualLoc(int i, int j)
        {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return fabs(eqn.A * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3 + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            //int si, sj;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                    {
                        resi = temp;
                        //si = i;
                        //sj = j;
                    }
                }
            //cout << "!!!!!!!!!!!!" << si << "," << sj << endl;
            return resi;
        }
        void execute();
    };

    template<>
    void gaussSeidelLineRelaxation<D_X>::execute()
    {
        int n_mtr = ie - is + 2;
        VectorXd x(n_mtr), matB(n_mtr);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        int i;
        for (int j = js;j < je;j++)
        {
            tripletList.clear();
            i = 0;
            if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(i, i, 1.0));
                tripletList.push_back(T(i, i + 1, -1.0));
                /*matA(i, i) = 1.0;
                matA(i, i + 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(i) = -(m.x()[i + 1] - m.x()[i]) * np->f(m.x()[i], m.y()[j]);

            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(i, i, 1.0));
                /*
                 matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(i) = dp->f(m.x()[i], m.y()[j]);
            }
            for (i++;i < ie;i++)
            {
                double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
                double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
                double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
                double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
                double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
                double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
                tripletList.push_back(T(i, i, eqn.A * (DXI1 + DXI2) * DXI3 + (DYI1 + DYI2) * DYI3));
                tripletList.push_back(T(i, i - 1, -eqn.A * DXI2 * DXI3));
                tripletList.push_back(T(i, i + 1, -eqn.A * DXI1 * DXI3));
                /*
                matA(i, i) = eqn.A * (DXI1 + DXI2) * DXI2 + (DYI1 + DYI2) * DYI3;
                matA(i, i - 1) = -eqn.A * DXI2 * DXI2;
                matA(i, i + 1) = -eqn.A * DXI1 * DXI2;
                */
                matB(i) = (m.now()[i][j + 1] * DYI1 + m.now()[i][j - 1] * DYI2) * DYI3;

            }
            if (m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(i, i, 1.0));
                tripletList.push_back(T(i, i - 1, -1.0));
                /*
                matA(i, i) = 1.0;
                matA(i, i - 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(i) = (m.x()[i] - m.x()[i - 1]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(i, i, 1.0));
                /*
                matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(i) = dp->f(m.x()[i], m.y()[j]);
            }
            SparseMatrix<double> matA(n_mtr, n_mtr);
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
            lscg.compute(matA);
            x = lscg.solve(matB);
            */

            BiCGSTAB<SparseMatrix<double> > solver;
            solver.compute(matA);
            x = solver.solve(matB);


            //MatrixXd sol = matA.colPivHouseholderQr().solve(matB);
            /*
            SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
            solver.analyzePattern(matA);
            solver.factorize(matA);
            x = solver.solve(matB);
            */

            for (i = 0;i <= ie;i++)
                m.now()[i][j] = x(i);



        }
    }

    template<>
    void gaussSeidelLineRelaxation<D_Y>::execute()
    {
        int n_mtr = je - js + 2;
        VectorXd x(n_mtr), matB(n_mtr);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        int j;
        for (int i = is;i < ie;i++)
        {
            tripletList.clear();
            j = 0;
            if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j + 1, -1.0));
                /*matA(i, i) = 1.0;
                matA(i, i + 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = -(m.y()[j + 1] - m.y()[j]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                 matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            for (j++;j < je;j++)
            {
                double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
                double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
                double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
                double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
                double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
                double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
                tripletList.push_back(T(j, j, eqn.A * (DXI1 + DXI2) * DXI3 + (DYI1 + DYI2) * DYI3));
                tripletList.push_back(T(j, j - 1, -DYI2 * DYI3));
                tripletList.push_back(T(j, j + 1, -DYI1 * DYI3));
                /*
                matA(i, i) = eqn.A * (DXI1 + DXI2) * DXI2 + (DYI1 + DYI2) * DYI3;
                matA(i, i - 1) = -eqn.A * DXI2 * DXI2;
                matA(i, i + 1) = -eqn.A * DXI1 * DXI2;
                */
                matB(j) = eqn.A * (m.now()[i + 1][j] * DXI1 + m.now()[i - 1][j] * DXI2) * DXI3;

            }
            if (m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j - 1, -1.0));
                /*
                matA(i, i) = 1.0;
                matA(i, i - 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = (m.y()[i] - m.y()[i - 1]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            SparseMatrix<double> matA(n_mtr, n_mtr);
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
            lscg.compute(matA);
            x = lscg.solve(matB);
            */

            BiCGSTAB<SparseMatrix<double> > solver;
            solver.compute(matA);
            x = solver.solve(matB);

            //MatrixXd sol = matA.colPivHouseholderQr().solve(matB);
            /*
            SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
            solver.analyzePattern(matA);
            solver.factorize(matA);
            x = solver.solve(matB);
            */

            for (j = 0;j <= je;j++)
                m.now()[i][j] = x(j);
        }
    }


    class ADITwoSteps :public scheme
    {
    public:
        int scflag;
        ellipticEuqation& eqn;
        gaussSeidelLineRelaxation<D_Y> sc1;
        gaussSeidelLineRelaxation<D_X> sc2;
        ADITwoSteps(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, ellipticEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input), sc1(m_input, is_input, ie_input, js_input, je_input, eqn_input), sc2(m_input, is_input, ie_input, js_input, je_input, eqn_input), scflag(0) {}
        double residualLoc(int i, int j)
        {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            return fabs(eqn.A * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3 + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            //int si, sj;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                    {
                        resi = temp;
                        //si = i;
                        //sj = j;
                    }
                }
            //cout << "!!!!!!!!!!!!" << si << "," << sj << endl;
            return resi;
        }
        void execute()
        {
            if (scflag == 0)
                sc1.execute();
            else if (scflag == 1)
                sc2.execute();
            scflag = (scflag + 1) % 2;
        }
    };
}//elliptic

namespace transonicSmallDisturbance
{
    class gaussSeidelLineRelaxation :public scheme
    {
    public:
        transonicSmallDisturbanceEuqation& eqn;
        gaussSeidelLineRelaxation(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, transonicSmallDisturbanceEuqation& eqn_input) :scheme(m_input, is_input, ie_input, js_input, je_input), eqn(eqn_input) {}
        double residualLoc(int i, int j)
        {
            double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
            double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
            double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
            double DXI4 = 1 / (m.x()[i - 1] - m.x()[i - 2]);
            double DXI5 = 2 / (m.x()[i] - m.x()[i - 2]);
            double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
            double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
            double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
            double fun = 0;
            if (eqn.A(i, j) > 0)
                fun += eqn.A(i, j) * ((m.now()[i + 1][j] - m.now()[i][j]) * DXI1 - (m.now()[i][j] - m.now()[i - 1][j]) * DXI2) * DXI3;
            if (i > 1 && eqn.A(i - 1, j) < 0)
                fun += eqn.A(i - 1, j) * ((m.now()[i][j] - m.now()[i - 1][j]) * DXI2 - (m.now()[i - 1][j] - m.now()[i - 2][j]) * DXI4) * DXI5;

            return fabs(fun + ((m.now()[i][j + 1] - m.now()[i][j]) * DYI1 - (m.now()[i][j] - m.now()[i][j - 1]) * DYI2) * DYI3);
        }
        double residual()
        {
            double resi = 0;
            double temp;
            //int si, sj;
            for (int j = js;j < je;j++)
                for (int i = is;i < ie;i++)
                {
                    temp = residualLoc(i, j);
                    if (temp > resi)
                    {
                        resi = temp;
                        //si = i;
                        //sj = j;
                    }
                }
            //cout << "!!!!!!!!!!!!" << si << "," << sj << endl;
            return resi;
        }
        void execute();
    };
    void gaussSeidelLineRelaxation::execute()
    {
        int n_mtr = je - js + 2;
        VectorXd x(n_mtr), matB(n_mtr);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        int j;
        for (int i = is;i < ie;i++)
        {
            tripletList.clear();
            j = 0;
            if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j + 1, -1.0));
                /*matA(i, i) = 1.0;
                matA(i, i + 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = -(m.y()[j + 1] - m.y()[j]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                 matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            for (j++;j < je;j++)
            {
                double DXI1 = 1 / (m.x()[i + 1] - m.x()[i]);
                double DXI2 = 1 / (m.x()[i] - m.x()[i - 1]);
                double DXI3 = 2 / (m.x()[i + 1] - m.x()[i - 1]);
                double DXI4 = 1 / (m.x()[i - 1] - m.x()[i - 2]);
                double DXI5 = 2 / (m.x()[i] - m.x()[i - 2]);
                double DYI1 = 1 / (m.y()[j + 1] - m.y()[j]);
                double DYI2 = 1 / (m.y()[j] - m.y()[j - 1]);
                double DYI3 = 2 / (m.y()[j + 1] - m.y()[j - 1]);
                double a = 0, b = 0;
                if (eqn.A(i, j) > 0)
                {
                    a += eqn.A(i, j) * (DXI1 + DXI2) * DXI3;
                    b += eqn.A(i, j) * (m.now()[i + 1][j] * DXI1 + m.now()[i - 1][j] * DXI2) * DXI3;
                }
                if (i > 1 && eqn.A(i - 1, j) < 0)
                {
                    a += -eqn.A(i - 1, j) * DXI2 * DXI5;
                    b += eqn.A(i - 1, j) * (-m.now()[i - 1][j] * DXI2 - (m.now()[i - 1][j] - m.now()[i - 2][j]) * DXI4) * DXI5;
                }
                tripletList.push_back(T(j, j, a + (DYI1 + DYI2) * DYI3));
                tripletList.push_back(T(j, j - 1, -DYI2 * DYI3));
                tripletList.push_back(T(j, j + 1, -DYI1 * DYI3));
                /*
                matA(i, i) = eqn.A * (DXI1 + DXI2) * DXI2 + (DYI1 + DYI2) * DYI3;
                matA(i, i - 1) = -eqn.A * DXI2 * DXI2;
                matA(i, i + 1) = -eqn.A * DXI1 * DXI2;
                */
                matB(j) = b;

            }
            if (m.boudarymap()[i][j]->type() == isNeumann)
            {
                tripletList.push_back(T(j, j, 1.0));
                tripletList.push_back(T(j, j - 1, -1.0));
                /*
                matA(i, i) = 1.0;
                matA(i, i - 1) = -1.0;
                */
                Neumann* np = dynamic_cast<Neumann*>(m.boudarymap()[i][j]);
                matB(j) = (m.y()[i] - m.y()[i - 1]) * np->f(m.x()[i], m.y()[j]);
            }
            else if (m.boudarymap()[i][j] != nullptr && m.boudarymap()[i][j]->type() == isDirichlet)
            {
                tripletList.push_back(T(j, j, 1.0));
                /*
                matA(i, i) = 1.0;
                */
                Dirichlet* dp = dynamic_cast<Dirichlet*>(m.boudarymap()[i][j]);
                matB(j) = dp->f(m.x()[i], m.y()[j]);
            }
            SparseMatrix<double> matA(n_mtr, n_mtr);
            matA.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
            lscg.compute(matA);
            x = lscg.solve(matB);
            */

            BiCGSTAB<SparseMatrix<double> > solver;
            solver.compute(matA);
            x = solver.solve(matB);

            //MatrixXd sol = matA.colPivHouseholderQr().solve(matB);
            /*
            SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
            solver.analyzePattern(matA);
            solver.factorize(matA);
            x = solver.solve(matB);
            */

            for (j = 0;j <= je;j++)
                m.now()[i][j] = x(j);
        }
    }

}
#endif