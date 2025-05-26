#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

# define alpha 1.0
# define m 1.0
# define delta 0.1
# define N 50

double nn[7] = {0, 0.9, 1.0, 1.1, 1.5, 2.0, 5.0};
int iter = 0;

void pochodne(double, double*, double*);

double E(double*);

double WP(int, int);

void rk4_vec(double, double, int, double *,
    void (*)( double , double * , double *) );

void zapisz(double, double *, std::ofstream * out);

int main (){

    double * s ;

    double n = 5000;
    double dt = 0.02;
    void (* f )( double , double * , double *); // wskaźnik do funkcji

    f = pochodne ; // wskaźnik do funkcji
    s =( double *) malloc ( 2*(N+1) * sizeof ( double )); // tablica rozwiązań
    
    for(int k=0; k<7; k++){
        std::string name = {"lab4_atomy"};
        name += ("_" + std::to_string(k) + ".txt");
        // sciezka folderu do zapisania
        std::string folder = "/Users/michau/Documents/fizyka-teoretyczna/lab4_lancuch_atomow/wyniki/";
        std::string path = folder + name; 
        std::ofstream outW(path);

        // warunki początkowe :
        double t=0;
        for(int i=0; i<2*(N+1); i++)
            s[i] = WP(i, k);
        zapisz(0, s, &outW);
        if(iter>0){
            double omegaN = 2*sqrt(alpha/m)*abs(sin(nn[iter]*M_PI/(2*N)));
            double t_max = 20*2*M_PI/omegaN;
            n = t_max/dt;
        }
        // symulacja w czasie :
        for (int i=0; i<n; i++){
            rk4_vec (t , dt , 2*(N+1) , s , f );
            t = t + dt ;
        // zapis wynikow do pliku
            zapisz(t, s, &outW);
        }
        outW.close();
        std::cout << "zapisane " << k << std::endl;
        iter++;
    }

    return 0;
}


void pochodne ( double t , double *s , double * k){

    k[0] = 0;
    k[N] = 0;
    k[N+1] = 0;
    k[2*N+1]= 0;
    for(int i=1; i<N; i++){
        k[i] = s[N+1+i];
        k[N+1+i] = alpha/m*(s[i-1]-2*s[i]+s[i+1]);
    }
    if(iter>0){
        double F=0.01;
        double omegaN = 2*sqrt(alpha/m)*abs(sin(nn[iter]*M_PI/(2*N)));
        k[N+1+1] = alpha/m*(s[N+1]-2*s[N+2]+s[N+3]) + F/m*sin(omegaN*t);
    }
}


double T(double *s){
    double sum = 0;
    for(int i=0; i<N+1; i++){
        sum += pow(s[N+1+i],2);
    }
    return 0.5*m*sum;
}

double U(double *s){
    double sum = 0;
    for(int i=1; i<N+1; i++){
        sum += pow(s[i-1]-s[i]+delta,2);
    }
    return 0.5*alpha*sum;
}

double WP(int i, int typ){
    if(typ == 0){
        if(i<N+1)
            return delta*i + delta/3*exp(-pow(delta*i-0.5*delta*N,2)/(2*pow(3*delta,2)));
        else
            return 0;
    }
    else
        if(i<N+1)
            return delta*i;
        else
            return 0;
}

void rk4_vec ( double t , double dt , int n , double *s ,
    void (* f )( double , double * , double *) ){
        # define M 102
        static double k1[M], k2 [M], k3[M], k4[M], w[M];
        int i ;
        for(i=0; i<n; i++)
            w[i] = s[i];
        f (t ,w , k1 );
        for(i=0; i<n; i++) 
            w[i]= s[i]+ dt/2 * k1[i];
        f (t + dt/2, w, k2);
        for(i=0; i<n; i++) 
            w[i] = s[i] + dt/2 * k2[i];
        f (t + dt/2, w, k3);
        for(i=0; i<n; i++) 
            w [i] = s[i] + dt * k3[i];
        f (t + dt, w, k4);
        for(i=0; i<n; i++) 
            s[i] = s[i] + dt/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

void zapisz(double t, double* s, std::ofstream * out){

    *out << t;
    for(int i=0; i<2*(N+1); i++){
        if(i<N+1)
            *out << " " << s[i]-i*delta;
        else
            *out << " " << s[i];
    }
    *out << " "  << T(s) << " " << U(s) << std::endl;
}


