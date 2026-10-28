// Steady-state temperature profiles (spin chain and spin ladder) from the same BSc thesis.

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;
double A(double C,double tau,double cp,double cm,double kp,double k)
{
    return sqrt((1/(C*tau))*(cp*cm/kp+pow(cp,2)/k));
}
double B(double C,double tau,double cp,double k,double kp)
{
    return (pow(cp,2))/(C*tau*k*kp);
}
double tempmag(double cp,double cm,double C,double tau,double kp,double D,double v,double k,double Qp,double x,double T0)
{
    return T0-Qp*B(C,tau,cp,k,kp)*x/pow(A(C,tau,cp,cm,kp,k),2);//+Qp*(B(C,tau,cp,k,kp)/pow(A(C,tau,cp,cm,kp,k),2)-1/kp)*sinh(A(C,tau,cp,cm,kp,k)*x)/(A(C,tau,cp,cm,kp,k)*cosh(A(C,tau,cp,cm,kp,k)/2));
}
double tempph(double Qp,double km,double kp,double x,double T0)
{
    return T0-(Qp/(km+kp))*x;
}



int main()
{
    double const cp=2.8*pow(10,6);
    double const cm=1.5*pow(10,5);
    double C=cp+cm;
    double const tau=4.0*pow(10,-4);
    double const kp=1.0;
    double const D=kp/cp;
    double km=45;
    double v=0.86;
    //double const v=sqrt(km/cm*tau);
    double k=1/(C*v*v*tau);
    double Qp=1;
    double T0=0;

    double x0=-0.5;
    double dx=0.0001;

    ofstream stat;
    stat.open ("static.txt");

    cout<<A(C,tau,cp,cm,kp,k)<<' '<<B(C,tau,cp,k,kp)<<endl;


    for(int i{0};i<=10000;i++)
        {
            double x=x0+i*dx;
            stat<<x0+i*dx<<' '<<tempmag(cp,cm,C,tau,kp,D,v,k,Qp,x,T0)<<' '<<tempph(Qp,km,kp,x,T0)<<endl;
        }



    return 0;
}
