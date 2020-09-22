#ifndef EQUATION_H
#define EQUATION_H

class ellipticEuqation
{
    public:
    //A\phi_{xx}+ \phi_{yy}=0
    double A;
    ellipticEuqation(double A_input):A(A_input){}
};

class transonicSmallDisturbanceEuqation
{
    public:
    //A\phi_{xx}+ \phi_{yy}=0
    meshBase& msh;
    double A(int i,int j)
    {
        double phix=(msh.old()[i+1][j]-msh.old()[i-1][j])/(msh.x()[i+1]-msh.x()[i-1]);
        return 1-m_*m_-(1.4+1)*m_*m_*phix/v_;
    }
    double m_,v_;
    transonicSmallDisturbanceEuqation(double m_input,double v_input,meshBase& mesh_input):m_(m_input),v_(v_input),msh(mesh_input){}
};
#endif