#ifndef BASEMESH_H
#define BASEMESH_H
#include <vector>
#include "solverManger.H"
template<int mesh_num>
struct mesh :public meshBase
{
    int nx, ny;
    int newflag;
    int oldflag;
    solverManger* psm;
    void setSolverManger(solverManger* sm_input) { psm = sm_input; }
    solverManger* solMng() { return psm; };
    double (*f)(double, double);
    std::vector<double> x_, y_;
    std::vector<std::vector<double> > value_[mesh_num];
    std::vector<std::vector<boundary*> > boudarymap_;
    mesh(int input_x, int input_y) :nx(input_x), ny(input_y), x_(input_x), y_(input_y), newflag(0), oldflag(mesh_num - 1), boudarymap_(input_x) {
        for (int j = 0;j < mesh_num;j++)
        {
            value_[j].resize(input_x);
            for (int i = 0; i < input_x; i++)
                value_[j][i].resize(input_y);
        }
        for (int i = 0; i < input_x; i++)
            boudarymap_[i].resize(input_y);
    }
    void setinitf(double (*f_input)(double, double)) { f = f_input; }
    void init()
    {
        for (auto& ib : psm->bl)
            for (int i = ib->is;i < ib->ie;i++)
                for (int j = ib->js;j < ib->je;j++)
                    boudarymap_[i][j] = ib;

        for (int k = 0;k < mesh_num;k++)
            for (int i = 0;i < nx;i++)
                for (int j = 0;j < ny;j++)
                {
                    value_[k][i][j] = f(x_[i], y_[i]);
                }

    }
    std::vector<std::vector<double> >& now() { return value_[newflag]; }
    std::vector<std::vector<double> >& old() { return value_[oldflag]; }
    void swap() { int t;t = newflag;newflag = oldflag;oldflag = t; }
    std::vector<double>& x() { return x_; }
    std::vector<double>& y() { return y_; }
    std::vector<std::vector<boundary*> >& boudarymap() { return boudarymap_; };
    int nLayer() { return mesh_num; }
    std::vector<std::vector<double> >* value() { return value_; }
};

template<>
void mesh<1>::swap() {}
#endif