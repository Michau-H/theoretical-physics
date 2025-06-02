#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>


const int n = 30, m = 150, j_1 = 60, j_2 = 90,  IT_MAX = 5000;
const double deltaZ =0.1, deltaRo = 0.1, V0 = 10;

// double ro(int, int);

// double roP(int, int, double**);

void warBrzeg(double**);

double relaksacja(int, int, double**);

double S(double**);

void zapiszS(int, double, std::ofstream *);

void zapiszU(double**, std::ofstream *);

int main(){


    double** u = new double*[n+1];
    for (int i = 0; i < n+1; i++) {
        u[i] = new double[m+1]();
    }
    
    // double** ro_m = new double*[2*N+1];
    // for (int i = 0; i < 2*N+1; i++) {
    //     ro_m[i] = new double[2*N+1]();
    // }

    // double** d_m = new double*[2*N+1];
    // for (int i = 0; i < 2*N+1; i++) {
    //     d_m[i] = new double[2*N+1]();
    // }

    // sciezka folderu do zapisania
    std::string folder = "/Users/michau/Documents/fizyka-teoretyczna/lab5_rownanie_Poissona/wyniki/";

    std::string path500 = folder + "u_at_5000.txt"; 
    std::ofstream outW500(path500);

    std::string pathS = folder + "wart_S.txt"; 
    std::ofstream outWS(pathS);

    
    for(int iter=0; iter<IT_MAX; iter++){
        for(int i=1; i<n; i++){
            for(int j=1; j<m; j++){
                u[i][j] = relaksacja(i,j,u);
            } 
        }

        warBrzeg(u);
        
        zapiszS(iter, S(u), &outWS);
        std::cout<< S(u) << std:: endl;
    }

    zapiszU(u, &outW500);
    outW500.close();

    outWS.close();


    for (int i = 0; i < n+1; i++) {
        delete[] u[i];
        // delete[] ro_m[i];
        // delete[] d_m[i];
    }
    delete[] u;
    // delete[] ro_m;
    // delete[] d_m;
    

    return 0;
}

// double ro(int x, int y){
//     return exp(-(pow((x-x0-N),2)+pow(y-N, 2))/(d*d)) - exp(-(pow((x+x0-N),2)+pow(y-N, 2))/(d*d));
// }

// double roP(int i, int j, double** u){
//     return -(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] -4*u[i][j])/pow(dx,2);
// }

void warBrzeg(double** u){
    for(int j=0; j<j_1; j++)
        u[n][j] = V0;
    for(int j=j_1; j<j_2; j++)
        u[n][j] = 0;
    for(int j=j_2; j<m+1; j++)
        u[n][j] = V0;
    for(int i=1; i<n; i++){
        u[i][m] = u[i][m-1];
        u[i][0] = u[i][1];
    }
    for(int j=0; j<m+1; j++)
        u[0][j] = u[1][j];
}

double relaksacja(int i, int j, double** u){
    return  1.0/(2.0/pow(deltaRo,2) + 2.0/pow(deltaZ,2))
            * ((u[i+1][j] + u[i-1][j])/(pow(deltaRo,2)) 
            + (u[i+1][j] - u[i-1][j])/(i*2*pow(deltaRo,2))
            + (u[i][j+1] + u[i][j-1])/(pow(deltaZ,2)));
}

double S(double ** u){
    double suma=0;
    for(int i=1; i<n; i++){
        for(int j=1; j<m; j++){
            suma += pow((u[i+1][j] - u[i-1][j])/(2*deltaRo),2)
            + pow((u[i][j+1] - u[i][j-1])/(2*deltaZ),2);
        }
        suma *= deltaRo*i;
    }
    return M_PI*deltaRo*deltaZ*suma;
}


void zapiszS(int i, double s, std::ofstream * out){
    *out << i << " " << s << std::endl;
}

void zapiszU(double** u, std::ofstream * out){
    for(int i=0; i<n+1; i++){
        for(int j=0; j<m+1; j++){
            *out << u[i][j] << " ";
        }
        *out << std::endl;
    }
}


