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

namespace p5 {
    const double gamma = 1.4;
    const double rho_inf = 1;
    const double p_inf = 1;
    const double a_inf = sqrt(gamma * p_inf / rho_inf);
    const double m_inf = 0.908;
    const double v_inf = a_inf * m_inf;
    double dphidy(double x, double y)
    {
        if (x < 0 || x>1)
            return 0;
        return -v_inf * 300 * (2 * x - 1) / sqrt(6205081 + 360000 * x - 360000 * x * x);
    }
    double boundryPhi(double x, double y)
    {
        return x * v_inf;
    }
}

namespace p6 {
    const double gamma = 1.4;
    const double rho_inf = 1;
    const double m_inf = 0.735;
    const double v_inf = 1.0;
    const double p_inf = 1.0 / (gamma * m_inf * m_inf);
    const double a_inf = v_inf / m_inf;

    double dphidy(double x, double y)
    {
        if (x < 0 || x>1)
            return 0;
        return -v_inf * 300 * (2 * x - 1) / sqrt(6205081 + 360000 * x - 360000 * x * x);
    }
    double boundryPhi(double x, double y)
    {
        return x * v_inf;
    }
}
double zero(double, double)
{
    return 0.0;
}
/*
double boundryPhi(double x, double y)
{
    return x * v_inf;
}
*/
/*
double dphidy(double x, double y)
{
    double a_inf = sqrt(1.4);
    double v_inf = a_inf * 0.5;
    if (x < 0 || x>1)
        return 0;
    return -v_inf * 300 * (2 * x - 1) / sqrt(6205081 + 360000 * x - 360000 * x * x);
}*/
int main()
{
    int gx1 = 15, gx2 = 20, gx3 = 15;
    int gy1 = 50;

    int tolgx = gx1 + gx2 + gx3;
    
    mesh<1>  m(tolgx + 1, gy1 + 1);

    m.x()[0] = -50;
    m.x()[gx1] = 0;
    m.x()[gx1 + gx2] = 1;
    m.x()[tolgx] = 51;

    geometric(m.x().begin(), m.x().begin() + gx1 + 1, 1.0 / gx2);
    linear(m.x().begin() + gx1, m.x().begin() + gx1 + gx2 + 1);
    geometric(1.0 / gx2, m.x().begin() + gx1 + gx2, m.x().begin() + tolgx + 1);

    m.y()[0] = 0;
    m.y()[gy1] = 50;
    geometric(0.006, m.y().begin(), m.y().end());
    m.setinitf(zero);

    Dirichlet left(m, 0, 1, 0, gy1 + 1, zero);
    Dirichlet top(m, 0, tolgx + 1, 50, gy1 + 1, zero);
    Dirichlet right(m, tolgx, tolgx + 1, 0, gy1 + 1, zero);
    //Dirichlet bottom(m, 1, 50, 0, 1, boundryPhi);
    Neumann bottom(m, 1, tolgx, 0, 1, 0, 1, p6::dphidy);
    left.init();
    top.init();
    right.init();
    bottom.update();
    transonicSmallDisturbanceEuqation eqn(p6::m_inf, p6::v_inf, m);

    //elliptic::pointJacobi sc(m, 1, 50, 1, 50, eqn);
    //elliptic::pointGaussSeidel sc(m, 1, 50, 1, 50, eqn);
    //elliptic::lineJacobi<D_X> sc(m, 1, 50, 1, 50, eqn);
    //elliptic::gaussSeidelLineRelaxation<D_Y> sc(m, 1, 50, 1, 50, eqn);
    transonicSmallDisturbance::gaussSeidelLineRelaxation sc(m, 1, tolgx, 1, gy1, eqn);
    solverManger solver(m, sc);
    solver.setBoundary(left);
    solver.setBoundary(top);
    solver.setBoundary(right);
    solver.setBoundary(bottom);
    solver.setTol(1e-5);
    solver.solve();
/*
    cin.get();

    transonicSmallDisturbanceEuqation eqn2(p5::m_inf, p5::v_inf, m);
    transonicSmallDisturbance::gaussSeidelLineRelaxation sc2(m, 1, tolgx, 1, gy1, eqn2);
    solverManger solver2(m, sc2);
    solver2.setBoundary(left);
    solver2.setBoundary(top);
    solver2.setBoundary(right);
    Neumann bottom2(m, 1, tolgx, 0, 1, 0, 1, p5::dphidy);
    solver2.setBoundary(bottom2);
    solver2.setTol(1e-6);

    for (int i = 0;i <= tolgx;i++)
        for (int j = 0;j <= gy1;j++)
            m.now()[i][j] *= p5::v_inf;

    solver2.solve();
    */







    /*
    double re=sc.residual();
    for (int i = 0;re > 1e-4;i++)
    {
        re=sc.residual();
        cout << i << " " << re << endl;
        sc.execute();
        bottom.update();
    }*/

    vector<double> u(tolgx - 1), v(tolgx - 1), p(tolgx - 1), cp(tolgx - 1);
    for (int i = 0;i < tolgx - 1;i++)
    {
        u[i] = p6::v_inf + (m.now()[i + 2][0] - m.now()[i][0]) / (m.x()[i + 2] - m.x()[i]);
        v[i] = (m.now()[i + 1][1] - m.now()[i + 1][0]) / (m.y()[1] - m.x()[0]);
        p[i] = p6::p_inf * pow(1 - (p6::gamma - 1) / 2.0 * p6::m_inf * p6::m_inf * ((u[i] * u[i] + v[i] * v[i]) / (p6::v_inf * p6::v_inf) - 1), p6::gamma / (p6::gamma - 1));
        cp[i] = (p[i] - p6::p_inf) / (0.5 * p6::rho_inf * p6::v_inf * p6::v_inf);
        // cout << i << " " << m.y[i] << endl;
    }




    /*    for (int i = 0; i < 51; i++)
            for (int j = 0; j < 51; j++) {
                m.value[i][j] = i * 51 + j;
            }
    */
    /*
        for (int i = 0; i < 51; i++)
            for (int j = 0; j < 1; j++)
            {
                cout << m.value[i][j] << endl;
            }
            */
 
    ofstream fout("cp6.csv");
    for (int i = 0; i < tolgx - 1; i++)
    {
        fout << m.x()[i + 1] << "," << u[i] << "," << v[i] << "," << cp[i] << "," << m.now()[i][0] << endl;
    }



    cin.get();
    return 0;
}
