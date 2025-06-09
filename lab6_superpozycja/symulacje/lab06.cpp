#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

const double mi0=1, eps0=1, R=1, sigma=1, L=3;
const int N=201, M=201, K=42;
const double dTheta = M_PI/(N-1), dPhi = 2*M_PI/(M-1);
const double dx=2*L/(K-1), dy=2*L/(K-1), omega0=0.1;

const double C=sigma*R*R*dTheta*dPhi/(4*M_PI*3*3);

int wspolczynniki(int, int);

void il_wekt(double*, double *, double*);

double il_skal(double*, double*);

void zapisz(double *, std::ofstream * out);

void zapiszU(double**, std::ofstream *);

int main(){
    double alpha[2] = {0, M_PI/4};

    double a[N], b[M];
    for(int i=0; i<N; i++){
        a[i] = wspolczynniki(i, N);
        
        b[i] = wspolczynniki(i, M);
        // std::cout<<a[i]<<" "<<b[i]<<std::endl;
    }
    

    std::string folder = "/Users/michau/Documents/fizyka-teoretyczna/lab6_superpozycja/wyniki/";
    

    for(int ii=0; ii<2; ii++){
        std::string name[5] = {"V", "EX", "EY", "BX", "BY"};
        
        std::string path = folder + name[0] + "_" + std::to_string(ii) + ".txt"; 
        std::ofstream outV(path);

        path = folder + name[1] + "_" + std::to_string(ii) + ".txt"; 
        std::ofstream outEX(path);
        path = folder + name[2] + "_" + std::to_string(ii) + ".txt"; 
        std::ofstream outEY(path);

        path = folder + name[3] + "_" + std::to_string(ii) + ".txt"; 
        std::ofstream outBX(path);
        path = folder + name[4] + "_" + std::to_string(ii) + ".txt"; 
        std::ofstream outBY(path);

        double** v = new double*[K];
        for (int i = 0; i < K; i++) {
            v[i] = new double[K]();
        }
        
        double** ex = new double*[K];
        for (int i = 0; i < K; i++) {
            ex[i] = new double[K]();
        }
        double** ey = new double*[K];
        for (int i = 0; i < K; i++) {
            ey[i] = new double[K]();
        }

        double** bx = new double*[K];
        for (int i = 0; i < K; i++) {
            bx[i] = new double[K]();
        }
        double** by = new double*[K];
        for (int i = 0; i < K; i++) {
            by[i] = new double[K]();
        }

        double omega[3]={omega0*sin(alpha[ii]), omega0*cos(alpha[ii]), 0};
        double x,y;
        for(int xi=0; xi<K; xi++){
            x=xi*dx-L;
            for(int yj=0; yj<K; yj++){
                y=yj*dy-L;
                double z=0;
                // double v=0, ex=0, ey=0, bx=0, by=0;
                for(int i=0; i<N; i++){
                    for(int j=0; j<M; j++){
                        double theta = dTheta*i;
                        double phi = dPhi*j;
                        double Rp[3] = {R*sin(theta)*cos(phi), R*sin(theta)*sin(phi),R*cos(theta)};
                        double rr[3]={x-Rp[0], y-Rp[1], z-Rp[2]};
                        double rr_mod = sqrt(il_skal(rr,rr));

                        v[xi][yj] += C/eps0*a[i]*b[j]*sin(theta)/rr_mod;
                        // std::cout<<rr[0]<<" "<<std::endl;
                        ex[xi][yj] += C/eps0*a[i]*b[j]*sin(theta)*rr[0]/pow(rr_mod,3);
                        ey[xi][yj] += C/eps0*a[i]*b[j]*sin(theta)*rr[1]/pow(rr_mod,3);

                        double g1[3], g2[3];
                        il_wekt(omega, Rp, g1);
                        il_wekt(rr, g1, g2);
                        bx[xi][yj] -= C*mi0*a[i]*b[j]*sin(theta)*g2[0]/pow(rr_mod,3);
                        by[xi][yj] -= C*mi0*a[i]*b[j]*sin(theta)*g2[1]/pow(rr_mod,3);
                    }
                }
                // double s[7]={x,y,v,ex,ey,bx,by};
                // zapisz(s, &outW);
                
            }
        }

        zapiszU(v, &outV);
        zapiszU(ex, &outEX);
        zapiszU(ey, &outEY);
        zapiszU(bx, &outBX);
        zapiszU(by, &outBY);

        for (int i = 0; i < K; i++) {
            delete[] v[i];
            delete[] ex[i];
            delete[] ey[i];
            delete[] bx[i];
            delete[] by[i];
        }
        delete[] v;
        delete[] ex;
        delete[] ey;
        delete[] bx;
        delete[] by;

        std::cout << "zapisane " << ii << std::endl;
    }

    return 0;
}


int wspolczynniki (int i, int N){
    int k;
    if ( i ==0 || i == N ) k =1;
    else if ( i %2==1 ) k =4;
    else if ( i %2==0 ) k =2;
    return k;
}

void il_wekt(double* c, double * d, double* g){
    g[0] = c[1]*d[2] - c[2]*d[1];
    g[1] = c[2]*d[0] - c[0]*d[2];
    g[2] = c[0]*d[1] - c[1]*d[0];
}

double il_skal(double* c, double* d){
    return c[0]*d[0] + c[1]*d[1] + c[2]*d[2];
}

void zapisz(double* s, std::ofstream * out){

    for(int i=0; i<7; i++){
        *out << s[i] << " ";
    }
    *out << std::endl;
}

void zapiszU(double** u, std::ofstream * out){
    for(int i=0; i<K; i++){
        for(int j=0; j<K; j++){
            *out << u[i][j] << " ";
        }
        *out << std::endl;
    }
}
