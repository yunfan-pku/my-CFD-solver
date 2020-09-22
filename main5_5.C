#include <iostream>
#include <vector>
#include <cmath>
#include "linesplit.H"

#include "mesh.H"
#include <fstream>
#include "Eigen/Core"
using namespace std;
using namespace linesplit;
using namespace Eigen;

const double gamma = 1.4;
const double rho_inf = 1;
const double p_inf = 1;
const double a_inf = sqrt(gamma * p_inf / rho_inf);
const double m_inf = 0.5;
const double v_inf = a_inf * m_inf;

double zero(double, double)
{
    return 0;
}
double boundryPhi(double x, double y)
{
    return x * v_inf;
}
double dphidy(double x, double y)
{
    double a_inf = sqrt(1.4);
    double v_inf = a_inf * 0.5;
    if (x < 0 || x>1)
        return 0;
    return -v_inf * 300 * (2 * x - 1) / sqrt(6205081 + 360000 * x - 360000 * x * x);
}
double dphidy6(double x, double y)
{
    double a_inf = sqrt(1.4);
    double v_inf = a_inf * 0.5;
    if (x < 0 || x>1)
        return 0;
    return -v_inf * 75 * (2 * x - 1) * 0.5 / sqrt(23716 + 5625 * x - 5625 * x * x);
}
int main()
{
    int gx1 = 15, gx2 = 20, gx3 = 15;
    int tolgx = gx1 + gx2 + gx3;
    int gy1 = 50;
    mesh<1>  m(gx1 + gx2 + gx3 + 1, gy1 + 1);

    m.x()[0] = -50;
    m.x()[gx1] = 0;
    m.x()[gx1 + gx2] = 1;
    m.x()[gx1 + gx2 + gx3] = 51;

    geometric(m.x().begin(), m.x().begin() + gx1 + 1, 1.0 / gx2);
    linear(m.x().begin() + gx1, m.x().begin() + gx1 + gx2 + 1);
    geometric(1.0 / gx2, m.x().begin() + gx1 + gx2, m.x().begin() + gx1 + gx2 + gx3 + 1);

    m.y()[0] = 0;
    m.y()[gy1] = 50;
    geometric(0.006, m.y().begin(), m.y().end());
    m.setinitf(zero);//boundryPhi);

    Dirichlet left(m, 0, 1, 0, gy1 + 1, zero);//boundryPhi);
    Dirichlet top(m, 0, tolgx + 1, gy1, gy1 + 1, zero);//boundryPhi);
    Dirichlet right(m, tolgx, tolgx + 1, 0, gy1 + 1, zero);//boundryPhi);
    //Dirichlet bottom(m, 1, 50, 0, 1, boundryPhi);
    Neumann bottom(m, 1, tolgx, 0, 1, 0, 1, dphidy);
    left.init();
    top.init();
    right.init();
    bottom.update();
    ellipticEuqation eqn(1 - m_inf * m_inf);


    vector<scheme*> sct;
    sct.push_back(new elliptic::pointJacobi(m, 1, tolgx, 1, gy1, eqn));
    sct.push_back(new elliptic::pointGaussSeidel(m, 1, tolgx, 1, gy1, eqn));
    sct.push_back(new elliptic::lineJacobi<D_Y>(m, 1, tolgx, 1, gy1, eqn));
    sct.push_back(new elliptic::gaussSeidelLineRelaxation<D_Y>(m, 1, tolgx, 1, gy1, eqn));
    sct.push_back(new elliptic::ADITwoSteps(m, 1, tolgx, 1, gy1, eqn));

    //elliptic::pointJacobi sc(m, 1, 50, 1, 50, eqn);
    //elliptic::pointGaussSeidel sc(m, 1, 50, 1, 50, eqn);
    //elliptic::lineJacobi<D_X> sc(m, 1, 50, 1, 50, eqn);
    //elliptic::gaussSeidelLineRelaxation<D_Y> sc(m, 1, 50, 1, 50, eqn);
    //elliptic::ADITwoSteps sc(m, 1, tolgx, 1, gy1, eqn);

    solverManger solver(m);
    solver.setBoundary(left);
    solver.setBoundary(top);
    solver.setBoundary(right);
    solver.setBoundary(bottom);
    solver.setTol(1e-6);

    vector<vector<double>> toltable(5);
    for (int i = 0;i < 5;i++)
    {

        solver.setScheme(*sct[i]);
        solver.solve(toltable[i]);
    }

    vector<double> u(tolgx - 1), v(tolgx - 1), p(tolgx - 1), cp(tolgx - 1);
    for (int i = 0;i < tolgx - 1;i++)
    {
        u[i] = v_inf + (m.now()[i + 2][0] - m.now()[i][0]) / (m.x()[i + 2] - m.x()[i]);
        v[i] = (m.now()[i + 1][1] - m.now()[i + 1][0]) / (m.y()[1] - m.x()[0]);
        p[i] = p_inf * pow(1 - (gamma - 1) / 2.0 * m_inf * m_inf * ((u[i] * u[i] + v[i] * v[i]) / (v_inf * v_inf) - 1), gamma / (gamma - 1));
        cp[i] = (p[i] - p_inf) / (0.5 * rho_inf * v_inf * v_inf);
    }

    /*
        for (int i = 0; i < 51; i++)
            for (int j = 0; j < 1; j++)
            {
                cout << m.value[i][j] << endl;
            }
    */
    ofstream fout("cp5-1.csv");
    for (int i = 0; i < tolgx - 1; i++)
    {
        fout << m.x()[i + 1] << "," << u[i] << "," << v[i] << "," << cp[i] << "," << m.now()[i][0] << endl;
    }
    ofstream fout2("res5-1.csv");
    for (int i = 0; i < 1000; i++)
    {
        for (int j = 0; j < 5; j++) {
            fout2<<toltable[j][i]<<",";
        }
        fout2<<endl;
    }

    cin.get();
    return 0;
}
