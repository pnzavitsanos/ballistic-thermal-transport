// Ballistic magnetic heat transport coupled to diffusive phonons.
// Time evolution solver from the BSc thesis "Ballistic thermal transport coupled to phonons"
// (P. Zavitsanos, University of Crete, 2020), basis of Phys. Rev. B 104, 024429 (2021).
//
// Method: Fourier method of lines (cosine/sine series on [0,1]), Simpson integration for the
// initial Fourier coefficients, fourth-order Runge-Kutta in time, inverse transform for profiles.
// Output files: phtemp.txt, magtemp.txt, system.txt, energy.txt, temperature.txt, profile.txt
//
// Usage: ./coupled_transport [D v cp cm]   (defaults: 0.1 1.0 0.9 0.1)

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
using namespace std;

//simple diffusion  :Limit 1/tau-->inf
double fkd(double x,double dqt)
{
    return -dqt*x;
}
void ruku_nd(int n,double h,double pi,vector <double> & acd,double dqt,double k1,double k2,double k3,double k4,double d)
{

    for(int i{0};i<n;i++)
        {
            dqt=d*pow(pi*i,2);
            k1=h*fkd(acd[i],dqt);
            k2=h*fkd(acd[i]+k1/2,dqt);
            k3=h*fkd(acd[i]+k2/2,dqt);
            k4=h*fkd(acd[i]+k3,dqt);
            acd[i]=acd[i]+(1.0/6.0)*(k1+2*k2+2*k3+k4);
        }
        return;
}

//Runge Kutta
double fk(double x,double y,double dqt,double cm)
{
    return -dqt*x-cm*(x-y);
}
double fl(double x,double y,double z,double vqt,double cp)
{
    return -vqt*z-cp*(y-x);
}
double fm(double y,double z,double vqt,double cp)
{
    return vqt*y-cp*(z);
}
void ruku_n(int n,double h,double pi,vector <double> & ac,vector <double> & bc,
vector <double> & cc,double cm,double cp,double vqt,double dqt,
double k1,double k2,double k3,double k4,double l1,double l2,double l3,
double l4,double m1,double m2,double m3,double m4,double d,double v)
{

    for(int i{1};i<n;i++)
        {
            dqt=d*pow(pi*i,2);
            vqt=v*pi*i;


            k1=h*fk(ac[i],bc[i],dqt,cm);
            l1=h*fl(ac[i],bc[i],cc[i],vqt,cp);
            m1=h*fm(bc[i],cc[i],vqt,cp);

            k2=h*fk(ac[i]+k1/2,bc[i]+l1/2,dqt,cm);
            l2=h*fl(ac[i]+k1/2,bc[i]+l1/2,cc[i]+m1/2,vqt,cp);
            m2=h*fm(bc[i]+l1/2,cc[i]+m1/2,vqt,cp);

            k3=h*fk(ac[i]+k2/2,bc[i]+l2/2,dqt,cm);
            l3=h*fl(ac[i]+k2/2,bc[i]+l2/2,cc[i]+m2/2,vqt,cp);
            m3=h*fm(bc[i]+l2/2,cc[i]+m2/2,vqt,cp);

            k4=h*fk(ac[i]+k3,bc[i]+l3,dqt,cm);
            l4=h*fl(ac[i]+k3,bc[i]+l3,cc[i]+m3,vqt,cp);
            m4=h*fm(bc[i]+l3,cc[i]+m3,vqt,cp);

            ac[i]=ac[i]+(1.0/6.0)*(k1+2*k2+2*k3+k4);
            bc[i]=bc[i]+(1.0/6.0)*(l1+2*l2+2*l3+l4);
            cc[i]=cc[i]+(1.0/6.0)*(m1+2*m2+2*m3+m4);
        }
        return;
}

void time_evolution0(vector <double> & ac,vector <double> & bc,vector <double> & cc,
double a0,double b0,double c0,double t,double cp,double cm)
{
    ac[0]=a0-cm*(a0-b0)*(1.0-exp(-t));
    bc[0]=b0-cp*(b0-a0)*(1.0-exp(-t));
    cc[0]=c0;
}

//Simpson integration
double f(double x,int n,double pi)
{
    return 2*exp(-300.0*(x-0.5)*(x-0.5))*cos(pi*n*x);
}
double intergal(double a,double b,int divisions,int n,double pi)
{
    double h=(b-a)/divisions;
    double s=(f(a,n,pi)+f(b,n,pi));
    for(int i{1};i<=divisions-1;i+=2)
        {
            s+= 4*f(a+i*h,n,pi);
        }
    for(int i{2};i<=divisions-2;i+=2)
        {
            s+=2*f(a+i*h,n,pi);
        }

    return (h/3)*s;
}



int main(int argc, char* argv[])
{
    const double pi=3.141592653589793238463;
    int constexpr n=17;                                // number of fourier const.
    int xmax=100.0;
    double dx=0.01;
    double x0=0.0;
    int m=10000;  //time steps
    double h=0.001;
    double dqt{0},vqt{0};
    double d{0.1};
    double v{1.0};
    double cm=0.1;
    double cp=0.9;
    double k1{0},k2{0},k3{0},k4{0},l1{0},l2{0},l3{0},l4{0},m1{0},m2{0},m3{0},m4{0};
    double t=0;
    if(argc==5)
    {
        d=atof(argv[1]); v=atof(argv[2]); cp=atof(argv[3]); cm=atof(argv[4]);
    }

    vector <double> ac(n),bc(n,0),cc(n,0),acd(n);

    //simpson
    int divisions=100;
    double a=0.0;
    double b=1.0;
    for(int k{0};k<n;k+=1)
        {
            ac[k]=intergal(a,b,divisions,k,pi);
            acd[k]=intergal(a,b,divisions,k,pi);
        }

    //iversion
    double const a0=ac[0];
    double const b0=bc[0];
    double const c0=0.0;
    
    double difftemp{0};
    double phtemp{0};
    double mptemp{0};
    double mmtemp{0};
    double energy{0},phenergy{0},magenergy{0};

    ofstream ph;
    ph.open ("phtemp.txt");
    ofstream mg;
    mg.open ("magtemp.txt");
    ofstream sy;
    sy.open ("system.txt");
    ofstream en;
    en.open ("energy.txt");
    ofstream pic;
    pic.open ("profile.txt");
    ofstream tem;
    tem.open ("temperature.txt");


    for(int j{0};j<m;j++)
        {
            t=h*j;
            if(j!=0)
            {
                time_evolution0(ac,bc,cc,a0,b0,c0,t,cp,cm);
                ruku_n(n,h,pi,ac,bc,cc,cm,cp,vqt,dqt,k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,d,v);

                ruku_nd(n,h,pi,acd,dqt,k1,k2,k3,k4,d);
            }
            for(int k{0};k<=xmax;k++)
                {
                    for(int i{0};i<n;i++)
                        {
                            if(i==0)
                            {
                                phtemp += 0.5*ac[i];
                                mptemp += 0.5*bc[i];
                                mmtemp += 0.5*bc[i];
                                difftemp +=0.5*acd[i];
                            }else{
                                phtemp += ac[i]*cos(i*pi*(x0+k*dx));
                                mptemp += bc[i]*cos(i*pi*(x0+k*dx))+cc[i]*sin(i*pi*(x0+k*dx));
                                mmtemp += bc[i]*cos(i*pi*(x0+k*dx))-cc[i]*sin(i*pi*(x0+k*dx));
                                difftemp +=acd[i]*cos(i*pi*(x0+k*dx));
                            }

                        }

                    ph<<x0+k*dx<<' '<<phtemp<<' '<<(difftemp-phtemp)<<' '<<t<<endl;
                    mg<<x0+k*dx<<' '<<0.5*(mptemp+mmtemp)<<' '<<t<<endl;
                    sy<<x0+k*dx<<' '<<phtemp+0.5*(mptemp+mmtemp)<<' '<<t<<endl;
                    
                    if(j>40){
                    if((j%50==0)&&(j<m))
                    {
                 		pic<<t<<' '<<x0+k*dx<<' '<<phtemp+0.5*(mptemp+mmtemp)
                 		<<' '<<0.5*(mptemp+mmtemp)<<' '<<phtemp<<' '<<abs(difftemp-phtemp)<<endl;
                    }
                    }

                    phenergy+=dx*cp*phtemp;
                    magenergy+=dx*cm*(mptemp+mmtemp)*0.5;
                    energy+=dx*(cp*phtemp+cm*0.5*(mptemp+mmtemp));

                    phtemp=0;
                    mptemp=0;
                    mmtemp=0;
                    difftemp=0;

                }

                en<<t<<' '<<phenergy<<' '<<magenergy<<' '<<energy<<endl;
                tem<<t<<' '<<phenergy/cp<<' '<<magenergy/cm<<' '<<phenergy/cp+magenergy/cm<<endl;
                
                ph<<endl<<endl;
                mg<<endl<<endl;
                sy<<endl<<endl;
                
                if(j>40){
                if((j%50==0)&&(j<m))
                {
                        pic<<endl<<endl;
                }
                }
                
                phenergy=0;
                magenergy=0;
                energy=0;
                
        }
return 0;
}

