#include<iostream>
#include<iomanip>
#include<math.h>
#include<cstdlib>
#include<fstream>
#include<string>
using namespace std;

double** array2d(int m, int n);
void save(double* data, int timestep);

int main(){
  int nx = 101;
  double* h;
  h = new double[nx+1];
  double* h0;
  h0 = new double[nx+1];
  double* wet;
  wet = new double[nx+1];
  double* eta;
  eta = new double[nx+1];
  double* etan;
  etan = new double[nx+1];
  double* u;
  u = new double[nx+1];
  double* un;
  un = new double[nx+1];
  double* w;
  w = new double[nx+1];
  double* wu;
  wu = new double[nx+1];
  double dt,dx,g,eps;
  double hmax,time,dtmax;
  double Apadle,Tpadle,pi,c,lambda;
  double dist,dwdz;
  int n,ntot,nout,mode;

  // init parameters
  dx = 10;
  dt = 0.1;
  g = 9.81;

  for(int k=1; k<nx+1;k++){
    h0[k] = 10.0;
    eta[k] = 0.0;
    etan[k] = 0.0;
  }
  h0[0]=0;
  h0[nx+1]=0;
  for(int k=0; k<nx+2;k++){
    h[k] = h0[k]+eta[k];
    wet[k] = 1;
    if(h[k] <= 0){wet[k] = 0;}
    u[k] = 0.0;
    un[k] = 0.0;
  }

  eps = 0;

  //wave padle
  Apadle = 1.0;
  Tpadle = 20.0;
  pi = 4.*atan(1.);

  // mode = 2 dam-break scenario
  mode=2;
  if(mode == 2){
    for(int k=46;k<57;k++){
      eta[k] = 1;
    }
  }

  ntot = 200.0/dt;
  nout = 2.0/dt;
  
  hmax = 0;
  for(int k=0; k<nx+2;k++){
    hmax = max(hmax,h[k]);
  }
  //max phase speed
  c = sqrt(g*hmax);
  lambda = dt*c/dx;
  if(lambda > 1){
    cout << "stability problem" << endl;
  }
  
  for(int n=0; n <ntot+1;n++){
    time=((double)n)*dt;
    if(mode==1){
      //wave padle at k = 51
      eta[51] = Apadle*sin(2.0*pi*time/Tpadle);
    }
    
   //dyn
    double pgradx; 
    double hue,huw,hwp,hwn,hen,hep;
    for(int k=0;k<nx+2;k++){
      if(wet[k]==1){
        pgradx = -g*(eta[k+1]-eta[k])/dx;
        un[k]=u[k]+dt*pgradx;
      }else{
        un[k]=0;
      }
    }
    for(int k=0;k<nx+2;k++){
      hep = 0.5*(un[k]+fabs(un[k]))*h[k];
      hen = 0.5*(un[k]-fabs(un[k]))*h[k+1];
      hue = hep+hen;
      hwp = 0.5*(un[k-1]+fabs(un[k-1]))*h[k-1];
      hwn = 0.5*(un[k-1]-fabs(un[k-1]))*h[k];
      huw = hwp+hwn;
      etan[k] = eta[k]-dt*(hue-huw)/dx;
    }
    //saphiro filter
    double term1,term2;
    for(int k=0;k<nx+2;k++){
      if(wet[k]==1){
        term1 = (1.0-0.5*eps*(wet[k+1]+wet[k-1]))*etan[k];
        term2 = 0.5*eps*(wet[k+1]*etan[k+1]+wet[k-1]*etan[k-1]);
        eta[k] = term1+term2;
      }else{
        eta[k] = etan[k];
      }
    }
    for(int k=0;k<nx+2;k++){
      h[k] = h0[k]+eta[k];
      u[k] = un[k];
    }

    // calculate w at 5 m below sea surface
    for(int k=1;k<nx+1;k++){
      dist = h[k]-5.0;
      w[k] = -(u[k]-u[k-1])/dx*dist;
    }
    w[0] = 0.;
    w[nx+1] = 0.;

    for(int k=1;k<nx+1;k++){
      wu[k] =0.5*(w[k]+w[k+1]);
    }
    if(n%100==0){
    cout << (double)n*dt << " " << eta[51]<< endl;
    save(eta,n); }
    
  }
}
void save(double* data, int timestep)
{
    ofstream myfile;
    char filename[11];
    sprintf(filename,"eta%d.dat", timestep);
    myfile.open(filename);
    std::locale system_locale("C");
    myfile.imbue(system_locale);

    for (int i = 0; i < 101+2; i++)
    {
            myfile << setprecision(10) << data[i] << "\n";
        }
    myfile.close();
    cout << filename << " saved!" << endl;
}
