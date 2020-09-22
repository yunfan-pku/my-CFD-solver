#ifndef BOUNDARY_H
#define BOUNDARY_H
#include"baseMesh.H"
#include<iostream>
using  namespace std;
const int isNeumann = 0;
const int isDirichlet = 1;
class boundary
{
public:
    meshBase& m;
    int is, ie;
    int js, je;
    boundary(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input) :m(m_input), is(is_input), ie(ie_input), js(js_input), je(je_input) {}
    virtual void init() = 0;
    virtual void update() = 0;
    virtual int type() = 0;
};

class Dirichlet : public boundary
{
public:
    int type() { return isDirichlet; }
    double (*f)(double, double);
    Dirichlet(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, double (*f_input)(double, double)) :boundary(m_input, is_input, ie_input, js_input, je_input), f(f_input) {}
    void update() {}
    void init() {
        for (int k = 0;k < m.nLayer();k++)
            for (int i = is;i < ie;i++)
                for (int j = js;j < je;j++)
                {
                    m.value()[k][i][j] = f(m.x()[i], m.y()[i]);
                }
    }
};
class Neumann : public boundary
{
public:
    int type() { return isNeumann; }
    double (*f)(double, double);
    int di, dj;
    Neumann(meshBase& m_input, int is_input, int ie_input, int js_input, int je_input, int di_input, int dj_input, double (*f_input)(double, double)) :boundary(m_input, is_input, ie_input, js_input, je_input), di(di_input), dj(dj_input), f(f_input) {}
    void update() {
        for (int i = is;i < ie;i++)
            for (int j = js;j < je;j++)
            {
                m.now()[i][j] = m.now()[i + di][j + dj] - f(m.x()[i], m.y()[i]) * (m.y()[j + dj] - m.y()[j] + m.x()[i + di] - m.x()[i]);

            }

    }
    void init() {}
};

#endif