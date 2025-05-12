#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>

# define g 9.81
# define R 1.0
# define m 1.0


void pochodne(double, double*, double*);

void rk4_vec(double, double, int, double *,
    void (*)( double , double * , double *) );

void zapisz(double, double *, std::ofstream * out);

int main (){
    double war[5]={4,45,90,135,175};
    double dt , * s ;
    int n , N ;
    void (* f )( double , double * , double *); // wskaźnik do funkcji

    n =2; dt =0.01;
    N = 1000; 

    f = pochodne ; // wskaźnik do funkcji
    s =( double *) malloc ( n * sizeof ( double )); // tablica rozwiązań
    
    for(int k=0; k<5; k++){
        std::string name = {"lab1_wahadlo"};
        name += ("_" + std::to_string(k) + ".txt");
        // sciezka folderu do zapisania
        std::string folder = "/Users/michau/Documents/fizyka-teoretyczna/lab1_wahadlo/wyniki/";
        std::string path = folder + name; 
        std::ofstream outW(path);
        // warunki początkowe :
        double t=0;
        s [0]=war[k]*M_PI/180;
        s [1]=0;
        zapisz(0, s, &outW);
        // symulacja w czasie :
        for (int i=0; i<N; i++){
            rk4_vec (t , dt , n , s , f );
            t = t + dt ;
        // zapis wynikow do pliku
            zapisz(t, s, &outW);
        }
        outW.close();
        std::cout << "zapisane " << k << std::endl;
    }

    return 0;
}


void pochodne ( double t , double *s , double * k ){
    k [0]= s [1];
    k [1]= - g / R * sin ( s [0]);
}



void rk4_vec ( double t , double dt , int n , double *s ,
    void (* f )( double , double * , double *) ){
    # define M 2
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
    *out << t << " " << s[0] << " " << s[1] << std::endl;
}


