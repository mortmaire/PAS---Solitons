#include <iostream>
#include <fftw3.h>
#include <cstdlib>
#include <complex>
#define complex complex<double>
#include <algorithm>
#include <stdint.h>
#include <iomanip>
// #include "gnuplot_i.hpp"
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
gest(gest),N(N){}

}p0;

class RANDOM{
    complex *tab;
    int n;
    bool is_cooked;
    int i;
public:
    RANDOM(int N):n(pow(2,N)),is_cooked(0),i(0){
        if(n==0)return;
        tab=new complex[n];
    }
private:
    double rand1(){
        double a = rand()/(double)RAND_MAX;
        if(a!=0)return a;
        else return 1;
 }
public:
    complex rand_n(){
        if(is_cooked==0){
            double r=sqrt(-log(rand1()));
            double t=2*M_PI*(rand1()+1/4.);
            tab[i]=polar(r,t);
        }
        if(i>=n-1){is_cooked=1;i=0;
//             random_shuffle(tab,tab+n-1);
        }
        return tab[i++];
    }
    void reset(){i=0;}
    void prepare(){
     while(is_cooked==0)rand_n();
     reset();
    }

}R(27);
class ffal{
public:
 int NZ;
 double *kk;
 complex* tab;
 parameters p;
 fftw_plan p1,p2;
 fftw_complex *in,*out;
 int TIME;
 double **zapis;
 int N;
 double gestosc;
 ffal(parameters pp=p0):p(pp){
    kk=new double[NZ];
//     zapis=new double[NZ];
    N=pow(2,ceil(log2(p.L/p.dx)));
    p.L=N*p.dx;
    tab=new complex[N];
    zapis = new double*[p.time];
    for(int t=0;t<p.time;t++)zapis[t]=new double[N];
//     zapis=new double[N*p.time];
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    p1  = fftw_plan_dft_1d(N,in,out, 1,FFTW_ESTIMATE);
    p2  = fftw_plan_dft_1d(N,in,out,-1,FFTW_ESTIMATE);
    TIME=0;
    R.reset();
 }
 ~ffal(){
  delete[] tab;
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_free(in);
  fftw_free(out);
//   for(int i=0;i<NZ;i++)delete[] zapis[i];
  delete[] kk;
  for(int t=0;t<p.time;t++)delete[] zapis[t];
  delete[] zapis;
 }
 inline complex eta(){return R.rand_n();}
 void evolve(){
    complex A(1,-p.gamma);
    complex A1=A*p.g*p.dt*complex(0,-1);
    complex A2=A*p.mu*p.dt*complex(0,1);
    complex *kin;
    double delta=sqrt(p.gamma*p.T*p.dt/p.dx);
    kin=new complex[N];
    for(int j=0;j<N;j++){
        double k;
        if(j<N/2)k=2*j*M_PI/p.L;
        else k=(j-N)*2*M_PI/p.L;
        kin[j]=exp(A*complex(0,-1)*k*k*p.dt)/(double)N;
    }

    for(int t=0;t<p.time;t++){
            zegar(t);
        for(int j=0;j<50;j++){
            for(int i=0;i<N;i++){
                static double f;
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
 }
 inline void zegar(int t){
  static int how_much=0;
  static double sr=0;
  static int ttime=time(0);
  how_much++;
  if(time(0)!=ttime){sr=how_much/(time(0)-ttime);ttime=time(0);how_much=0;}
  if(gpid==getpid())cerr<<"\r"<<p.T<<"\t"<<t<<"\t"<<sr<<flush;
 }
 
//  void fft(double **in0, double **out0,int w=-1){
//      for(int t=0;t<p.time;t++){
//          for(int i=0;i<N;i++){in[i][0]=in0[t][i];in[i][1]=in0[t][i];}
//          if(w==-1)fftw_execute(p2);
//          if(w==1)fftw_execute(p1);
//          for(int i=0;i<N;i++)
//          out0[t][i]=abs(complex(out[i][0],out[i][1]));
//      }
//  }
 
 void plot(double** buff=NULL,const char* opis="",int numb=0){
     if(buff==NULL)buff=zapis;
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
     sprintf(s,"T%04d.png",(int)p.T);
     image.write(s);
 }
 void bw(double **buff=NULL,double thr=0,double tt=0){

  if(buff==NULL)buff=zapis;
  if(tt!=0){
    for(int t=0;t<p.time;t++)for(int i=0;i<N;i++)if(buff[t][i]<tt)buff[t][i]=0;
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
 
 void save(const char* name="zapis.txt"){
  fstream file(name,ios::out);
  for(int i=0;i<N*p.time;i++)file<<zapis[i/N][i%N]<<endl;
  file.close();
 }
 void load(const char* name="zapis.txt"){
    fstream file(name,ios::in);
//     double z;
    for(int i=0;i<N*p.time;i++)
        file>>zapis[i/N][i%N];
    file.close();
    
        for(int t=p.time/2;t<p.time;t++)for(int i=0;i<N;i++)gestosc+=zapis[t][i];
    gestosc/=p.time*N/2;   
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
 
 
};


int main(){
    p0.time=2000;
    ffal psi;
    psi.evolve();
    psi.plot();
}



// int main(){
// 
//     p0.time=10000;
//     int N=p0.N;
//     fftw_plan p1,p2;
//     fftw_complex *in,*out;
//     in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
//     p1  = fftw_plan_dft_1d(N,out,in, 1,FFTW_ESTIMATE);
//     p2  = fftw_plan_dft_1d(N,in,out,-1,FFTW_ESTIMATE);
//     
// 
//     
//     
//     int n=10;
//     int x=1;
//     x=(fork()==0)*8+(fork()==0)*4+(fork()==0)*2+(fork()==0)+1;
//     if(x>n)return 0;
//     int gg=getpid();
//     
//     double **buff;
//     buff=new double*[p0.time];
//     for(int t=0;t<p0.time;t++)buff[t]=new double[N];
//     
//     double **buff2;
//     buff2=new double*[p0.time];
// 
//     
//     for(p0.T=x;p0.T<1000;p0.T+=n){
//         ffal psi;
//         psi.evolve();
//         char ss[64];
//         sprintf(ss,"../temp/T%04d.dat",(int)p0.T);
//         psi.save(ss);
// //             psi.load(ss);
//         for(int i=0;i<N*p0.time;i++)buff[i/N][i%N]=psi.zapis[i/N][i%N];
//         
//         sprintf(ss,"%d.txt",gg);
// //         fstream file(ss,ios::out|ios::app);
// //         file<<p0.T<<" "<<0<<" ";
// //         psi.plot(buff,"a",0);
// //         psi.bw(buff,0,psi.gestosc);
// //         psi.bw(buff,500);
// //         pair <double,double>a=psi.statystyka(buff);
// // 
// //         file<<a.first<<" "<<a.second<<" ";
// //         file.close();
//         
// //     for(double sigma=10;sigma<=500;sigma+=50){
// //         if(sigma==60)sigma=50;
// //         for(int tt=512;tt>0;tt/=2){
// //         for(int t=0;t<p0.time;t++){
// //          for(int i=0;i<N;i++){in[i][0]=psi.zapis[t][i];in[i][1]=0;}
// //          fftw_execute(p2);
// //          for(int i=tt;i<N/2;i++){          
// //              out[i][0]*=0;out[i][1]*=0;}
// //          for(int i=N/2;i<N-tt;i++){          
// //              out[i][0]*=0;out[i][1]*=0;}
// //          fftw_execute(p1);
// //         for(int i=0;i<N;i++)buff[t][i]=abs(complex(in[i][0],in[i][1]))/N;
// //      }
// //     }
//     for(double sigma=10;sigma<=500;sigma+=50){
//      if(sigma==60)sigma=50;
//         for(int t=0;t<p0.time;t++){
//          for(int i=0;i<N;i++){in[i][0]=psi.zapis[t][i];in[i][1]=0;}
//          fftw_execute(p2);
//          for(int i=0;i<N;i++){  
//             double m=(exp(-i*i/sigma/sigma/2)+exp(-(i-N)*(i-N)/sigma/sigma/2))/(1+exp(-N*N/sigma/sigma/2));            
//              out[i][0]*=m;
//              out[i][1]*=m;}
//          fftw_execute(p1);
//         for(int i=0;i<N;i++)buff[t][i]=abs(complex(in[i][0],in[i][1]))/N;
//         }
//         for(int proc=10;proc<200;proc+=10){
//         for(int t=0;t<p0.time;t++)
//         for(int i=0;i<N;i++)buff2[t][i]=buff[t][i];
//         sprintf(ss,"%d.txt",gg);
//         fstream file(ss,ios::out|ios::app);
//         psi.bw(buff2,0,psi.gestosc*proc/100.);
//         for(int i=100;i<=2000;i+=100){
//         psi.bw(buff2,i);
//         pair <double,double>a=psi.statystyka(buff2);
//         file<<p0.T<<" "<<proc<<" "<<sigma<<" "<<i<<" "<<a.first<<" "<<a.second<<"\n";
//         }
//         file.close();
//         }
//      }
//     }
// }
//  

    
//Kawałek z którego zrobiłem porównanie metody obcięcia do metody funkcji gaussa    
   /*
        int n=10;
        int x=(fork()==0)*8+(fork()==0)*4+(fork()==0)*2+(fork()==0)+1;
    p0.time=2000;
    int T0=1;
    if(x>n)return 0;
    for(p0.T=x;p0.T<=600;p0.T+=T0*n){
        
    if(p0.T>20&&T0<5){T0=5;p0.T=(p0.T-20)*5+20;}
    if(p0.T>50&&T0<10){T0=10;p0.T=(p0.T-50)*2+50;}
    if(p0.T>100&&T0<50){T0=50;p0.T=(p0.T-100)*5+100;}
  
//   cout<<"\\T{"<<setfill('0')<<setw(4)<<(int)p0.T<<"}{"<<setw(0)<<(int)p0.T<<"}"<<endl;
//     }
//   
    ffal psi;
    psi.evolve();
    
    
    
    char ss[64];
    sprintf(ss,"../temp/T%04d.dat",(int)p0.T);
    psi.save(ss);
//     psi.load(ss);
    double **buff;
    buff=new double*[p0.time];
    for(int i=0;i<p0.time;i++)buff[i]=new double[N];
    for(int t=0;t<p0.time;t++)for(int j=0;j<N;j++)buff[t][j]=psi.zapis[t][j];


    psi.plot(buff,"a",0);
    psi.bw(buff,0,psi.gestosc);
    psi.plot(buff,"b",0);
    psi.bw(buff,500);
    psi.plot(buff,"c",0);
    psi.bw(buff,2000);
    psi.plot(buff,"d",0);
        
    for(int i0=64;i0>0;i0-=4){

     for(int t=0;t<p0.time;t++){
         for(int i=0;i<N;i++){in[i][0]=psi.zapis[t][i];in[i][1]=0;}
         fftw_execute(p2);
         for(int i=i0;i<N-i0;i++){out[i][0]=0;out[i][1]=0;}
         fftw_execute(p1);
             for(int i=0;i<N;i++)buff[t][i]=abs(complex(in[i][0],in[i][1]))/N;
     }
    psi.plot(buff,"a",i0);
    psi.bw(buff,0,psi.gestosc);
    psi.plot(buff,"b",i0);
    psi.bw(buff,500);
    psi.plot(buff,"c",i0);
    psi.bw(buff,2000);
    psi.plot(buff,"d",i0);
    }
    double del=100;
    for(double sigma=500;sigma>=1;sigma-=del){
        if(sigma==100)del=10;
        if(sigma==20)del=5;
        if(sigma==10)del=1;
        
        
        for(int t=0;t<p0.time;t++){
         for(int i=0;i<N;i++){in[i][0]=psi.zapis[t][i];in[i][1]=0;}
         fftw_execute(p2);
         for(int i=0;i<N;i++){
             double m=(exp(-i*i/sigma/sigma/2)+exp(-(i-N)*(i-N)/sigma/sigma/2))/(1+exp(-N*N/sigma/sigma/2));            
             out[i][0]*=m;out[i][1]*=m;}
         fftw_execute(p1);
        for(int i=0;i<N;i++)buff[t][i]=abs(complex(in[i][0],in[i][1]))/N;
     }
    psi.plot(buff,"e",sigma);
    psi.bw(buff,0,psi.gestosc);
    psi.plot(buff,"f",sigma);
    psi.bw(buff,500);
    psi.plot(buff,"g",sigma);
    psi.bw(buff,2000);
    psi.plot(buff,"h",sigma);
        
        
    }
    
    }   }*/