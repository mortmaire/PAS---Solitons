#include <iostream>
#include <fftw3.h>
#include <cstdlib>
#include <complex>
#include <cmath>
#define complex complex<double>
#include <algorithm>
#include <iomanip>
#include <list>
#include <set>
#include <png++/png.hpp>

using namespace std;
using namespace png;

int seed=112233;
int gpid=getpid();

class parameters{
public:
 int T0;
 double dx;
 double L;
 double gamma;
 double T;
 double mu;
 double g;
 double dt;
 int time;
 double valmax;
 unsigned long int NB;
 double gest;
 int N;
 parameters(double dx=0.1,double L=200,double gamma=0.01,double T=1,
double mu=1,double g=0.01,double dt=0.01,double time=2000,double valmax=20,double gest=1.,double N=2048):
dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time),valmax(valmax),
gest(gest),N(N){
    N=pow(2,ceil(log2(L/dx)));
    L=N*dx;
}

}p0;



class ffal{
public:
    complex* tab;
    parameters p;
    fftw_plan p1,p2;
    fftw_complex *in,*out;
    int N;
    double gestosc;
    double **zapis;

    complex A,A1,A2;
    complex *kin;
    double delta;

    ffal(parameters pp=p0):p(pp){
        N=p.N;
        
        tab=new complex[N];
        zapis = new double*[p.time];
        for(int t=0;t<p.time;t++)zapis[t]=new double[N];
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
        p1  = fftw_plan_dft_1d(N,in,out, 1,FFTW_ESTIMATE);
        p2  = fftw_plan_dft_1d(N,in,out,-1,FFTW_ESTIMATE);
        
        A=complex(1,-p.gamma);
        A1=A*p.g*p.dt*complex(0,-1);
        A2=A*p.mu*p.dt*complex(0,1);
        
        delta=sqrt(p.gamma*p.T*p.dt/p.dx);
        kin=new complex[N];
        for(int j=0;j<N;j++){
            double k;
            if(j<N/2)k=2*j*M_PI/p.L;
            else k=(j-N)*2*M_PI/p.L;
            kin[j]=exp(A*complex(0,-1)*k*k*p.dt)/(double)N;
        }
        
    }
    ~ffal(){
        delete[] tab;
        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_free(in);
        fftw_free(out);
        for(int t=0;t<p.time;t++)delete[] zapis[t];
        delete[] zapis;
        delete[] kin;
    }
        inline double rand_1(){
        return (rand()+1.)/RAND_MAX;
 }
    
     inline complex eta(){
            double r=sqrt(-log(rand_1()));
            double t=2*M_PI*(rand_1()+1/4.);
            return polar(r,t);
         
    }
     
    inline void zegar(int t){
        static int how_much=0;
        static double sr=0;
        static int ttime=time(0);
        how_much++;
        if(time(0)!=ttime){sr=how_much/(time(0)-ttime);ttime=time(0);how_much=0;}
        if(gpid==getpid())cerr<<"\r"<<p.T<<"\t"<<t<<"\t"<<sr<<flush;
        }
     
    void evolve(int time=p0.time){
        if(time!=p.time){
            for(int t=0;t<p.time;t++)delete[] zapis[t];
            delete[] zapis;
            p.time=time;
            zapis=new double*[p.time];
            for(int t=0;t<p.time;t++)zapis[t]=new double[N];
        }

            for(int t=0;t<p.time;t++){
            zegar(t);
        for(int j=0;j<50;j++){
            for(int i=0;i<N;i++){
                static double f=1;
                f=norm(tab[i]);
                tab[i]*=exp(A1*f+A2);
                tab[i]+=delta*eta();
                f=norm(tab[i]);
                tab[i]*=exp(A1*f+A2);
                tab[i]+=delta*eta();
                }
            for(int i=0;i<N;i++){
                in[i][0]=tab[i].real();
                in[i][1]=tab[i].imag();
                }
            fftw_execute(p2);
            for(int i=0;i<N;i++){
                tab[i]=complex(out[i][0],out[i][1]);
                tab[i]*=kin[i];
                in[i][0]=tab[i].real();
                in[i][1]=tab[i].imag();
                }
            fftw_execute(p1);
            for(int i=0;i<N;i++)tab[i]=complex(out[i][0],out[i][1]);
            for(int i=0;i<N;i++)zapis[t][i]=abs(tab[i]);
            }
        }
        gestosc=0;
        for(int t=p.time/2;t<p.time;t++)for(int i=0;i<N;i++)gestosc+=zapis[t][i];
        gestosc/=p.time*N/2;
        cerr<<endl;
 }
 
 void plot(const char* opis=""){
     double** buff=zapis;
     image< rgb_pixel > image(p.time,N/2);
     double max=0;
     for(int t=0;t<p.time;t++)
         for(int i=0;i<N;i++)
             if(buff[t][i]>max)max=buff[t][i];
    bool is_bw=0;      
    double a[3];
    a[0]=buff[0][0];
    a[1]=a[0];
    for(int i=0;i<p.time*N;i++)
            if(buff[i/N][i%N]!=a[0]){a[1]=buff[i/N][i%N];break;}
    a[2]=a[1];
    for(int i=0;i<p.time*N;i++)
            if(buff[i/N][i%N]!=a[0] and buff[i/N][i%N]!=a[1])
            {a[2]=buff[i/N][i%N];break;}
    if(a[2]==a[1])is_bw=1;
    if(is_bw==0){
     for(int y=0;y<N/2;y++)
            for(int x=0;x<p.time;x++){
                double val=buff[x][2*y]/max;
                int r=sqrt(val)*255;
                int g=val*val*val*val*255;
                int b=val*(0.5-val)*16*255;
                if(b<0)b=0;
                image[y][x]=rgb_pixel(r,g,b);
            }
    }
    else{
        for(int y=0;y<N/2;y++)
            for(int x=0;x<p.time;x++){
                double val=buff[x][2*y]/max;
                image[y][x]=rgb_pixel(255*val,255*val,255*val);
            }
    }
     char s[64];
     sprintf(s,"T%04d%s.png",(int)p.T,opis);
     image.write(s);
 }
  void bw(double **buff=NULL,double thr=0,double tt=0){

  if(buff==NULL)buff=zapis;
  
  double aa1,aa2,aa3;
  aa1=0;
  aa2=aa1;
  aa3=aa1;
  for(int i=0;i<N*p.time;i++){
      if(aa1==aa2)if(buff[i/N][i%N]>aa1)aa2=buff[i/N][i%N];
      if(aa1!=aa2)if(buff[i/N][i%N]>aa2)aa3=buff[i/N][i%N];
//     if(aa1 != aa2 and aa2!=aa3)break;
  }
//   cout<<aa1<<" "<<aa2<<" "<<aa3<<endl;
  if(aa1!=aa2 and aa2!=aa3 and aa1!=aa3)
  if(tt!=0){
    for(int t=0;t<p.time;t++)for(int i=0;i<N;i++)
        if(buff[t][i]<tt)buff[t][i]=0;
        else buff[t][i]=1;
  }
  
  
    if(thr>0){
    for(int j=0;j<p.time;j++)
     for(int i=0;i<N;i++){
      int n=bucket(buff,i,j,0,2);
      if(n<thr)bucket(buff,i,j,2,1);
      else bucket(buff,i,j,2,3);
     }
    for(int j=0;j<p.time;j++)
        for(int i=0;i<N;i++)if(buff[j][i]==3)buff[j][i]=0;
    }
    
 }
 
 int bucket(double **buff,int i,int j,int from,int dest){
 set<pair<int,int> > zbior;
 if(buff[j][i]!=from)return 0;
 int n=0;
 zbior.insert(pair<int,int>(i,j));
 while(zbior.size()){
  set<pair<int,int> >::iterator it=zbior.begin();
  int x0=it->first;
  int y0=it->second;
  buff[y0][x0]=dest;
  zbior.erase(it);
  n++;
  for(int i=0;i<4;i++){
      int x;int y;
   if(i<2){x=x0+(i%2)*2-1;y=y0;}
   else{x=x0;y=y0+(i%2)*2-1;}
   x+=N;
   x%=N;
   if(y>= p.time or y<0)continue;
   if(buff[y][x]==from)zbior.insert(pair<int,int>(x,y));
  }
 }
  return n;
}
 
  pair<double,double> statystyka(double **buff=NULL){
        if(buff==NULL)buff=zapis;
        double n0=0;
        double d0=0;
        for(int t=p.time/2;t<p.time;t++){
            bool f=0;
            double n=0;
            double d=0;
                for(int i=0;i<N;i++){
                    int z=buff[t][i];
                    if(z==0 and f==0){d++;n++;f=1;}
                    else if(z==0 and f==1){d++;}
                    if(z==1 and f==1){f=0;}
            }
        d/=n;
        if(n==0)d=0;
        n0+=n;
        d0+=d;
    }
        n0/=p0.time/2;
        d0/=p0.time/2;
        return pair<double,double>(n0,d0);
 }
 pair<double,double> GetStatistics(double** buff=NULL,double part=1,int threshold=0){
     bool f=0;
     if(buff==NULL){
         f=1;
         buff=new double*[p.time];
         for(int t=0;t<p.time;t++){buff[t]=new double[N];
         for(int i=0;i<N;i++)buff[t][i]=zapis[t][i];
         }
    }
     bw(buff,threshold,gestosc*part);
     pair<double,double> ss=statystyka(buff);
     if(f){for(int i=0;i<p.time;i++)delete[] buff[i];
     delete[] buff;}
     return ss;
 }
};

int main(){
    ffal psi;
    psi.evolve(2000);
    for(int i=0;i<200;i++)
        for(int j=0;j<2000;j++){
            pair<double,double> ss=psi.GetStatistics(0,i/100.,j);
            cout<<i<<" "<<j<<" "<<ss.first<<" "<<ss.second<<endl;
        }
}