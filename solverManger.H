#ifndef SOLVERMANAGER_H
#define SOLVERMANAGER_H
#include "scheme.H"
#include<list>

class solverManger
{
public:
    meshBase& m;
    scheme* psc;
    std::list<boundary*> bl;
    double tol;
    solverManger(meshBase& m_input, scheme& sc_input) :m(m_input), psc(&sc_input), tol(1e-8) { psc->psm = this;m.setSolverManger(this); }
    solverManger(meshBase& m_input) :m(m_input), psc(nullptr), tol(1e-8) { m.setSolverManger(this); }
   
    void setScheme(scheme& sc_input) { psc=&sc_input;psc->psm = this; }
    void setBoundary(boundary& p) { bl.push_back(&p); }
    void setBoundary(boundary* p) { bl.push_back(p); }
    void setTol(double tol_input) { tol = tol_input; }
    void solve(vector<double>& table) {
        psc->init();
        m.init();
        for (auto& i : bl)
            i->init();
        for (auto& i : bl)
            i->update();
        double res = psc->residual();
        for (int i = 0;res > tol&&i<=1000;i++)
        {
            psc->execute();
            for (auto i : bl)
                i->update();
            res = psc->residual();
            table.push_back(res);
            if (i % 10 == 0)
                cout << i << " " << res << endl;
                
            m.swap();

        }
    }

};
#endif