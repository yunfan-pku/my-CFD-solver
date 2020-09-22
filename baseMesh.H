#ifndef MESHBASE_H
#define MESHBASE_H
#include <vector>
class boundary;
class solverManger;
struct meshBase
{
    virtual void setSolverManger(solverManger*) = 0;
    virtual solverManger* solMng() = 0;
    virtual std::vector<std::vector<boundary*> >& boudarymap() = 0;
    virtual std::vector<std::vector<double> >& now() = 0;
    virtual std::vector<std::vector<double> >& old() = 0;
    virtual void swap() = 0;
    virtual void init() = 0;
    virtual std::vector<double>& x() = 0;
    virtual std::vector<double>& y() = 0;
    virtual int nLayer() = 0;
    virtual std::vector<std::vector<double> >* value() = 0;
};


#endif

