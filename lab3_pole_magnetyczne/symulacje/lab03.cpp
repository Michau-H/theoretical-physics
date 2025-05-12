#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

# define q 1.0
# define m 1.0
# define B 1.0


void pochodne(double, double*, double*);

double E(double*);

void rk4_vec(double, double, int, double *,
    void (*)( double , double * , double *) );

void zapisz(double, double *, std::ofstream * out);

int main (){
    double war[2][4]={{1.1,1.0,0.0}, {M_PI/99,10,0,0}};
    double dt[2]={0.1,0.01};
    double * s ;
    int n =4;
    double N[2] ={500,5000};
    void (* f )( double , double * , double *); // wskaźnik do funkcji

    f = pochodne ; // wskaźnik do funkcji
    s =( double *) malloc ( n * sizeof ( double )); // tablica rozwiązań
    
    for(int k=0; k<2; k++){
        std::string name = {"lab2_stozek"};
        name += ("_" + std::to_string(k) + ".txt");
        // sciezka folderu do zapisania
        std::string folder = "/Users/michau/Documents/fizyka-teoretyczna/lab2_stozek/wyniki/";
        std::string path = folder + name; 
        std::ofstream outW(path);
        // warunki początkowe :
        double t=0;
        s [0]=war[k][0];
        s [1]=war[k][1];
        s [2]=war[k][2];
        s [3]=war[k][3];
        zapisz(0, s, &outW);
        // symulacja w czasie :
        for (int i=0; i<N[k]; i++){
            rk4_vec (t , dt[k] , n , s , f );
            t = t + dt[k] ;
        // zapis wynikow do pliku
            zapisz(t, s, &outW);
        }
        outW.close();
        std::cout << "zapisane " << k << std::endl;
    }

    return 0;
}


void pochodne ( double t , double *s , double * k ){

    k [0]= s[3]/m;
    k [1]= s[4]/(m*pow(s[0],2)) - q*B/(2*m);
    k [2]= s[5]/m;
    k [3]= pow(s[4],2)/(m*s[3])-pow(q*B,2)*s[0]/(4*m);
    k [4]= 0;
    k [5]= 0;
    return ;
    }

// double E(double *s){
//     return 0.5*(pow(tan(alfa)*s[1]*s[2], 2) + pow(s[3]/cos(alfa),2)) + g*s[1]*sin(alfa)*(1-cos(s[0]));
// }

void rk4_vec ( double t , double dt , int n , double *s ,
    void (* f )( double , double * , double *) ){
        # define M 6
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

    *out << t << " " << s[0] << " " << s[1] << " " << s[2] << " " << s[3] << " " << s[4] << " " << s[5] << std::endl;
}


