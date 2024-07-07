#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//n de cuerpos
#define C 7

/*
0 - sol
1 - terra
2 - mart
3 - jupiter
4 - mercurio
5- luna (tierra)
*/

//num de iteraciones
#define N 1051200

const double PI = 3.14159265359;
//LONGITUD PRÒPIA: 1U.A
const double L0 = 149597870.7; //km
//TEMPS PROPI: 1 Any/2Pi
const double T0 = (365.0*24.0*3600.0)/(2*PI);//s 

double MM[C]= { //masas
    1.99*pow(10.,30.),
    5.97*pow(10.,24.),
    6.42*pow(10.,23.),
    1.9*pow(10.,27.),
    3.302*pow(10,23),
    4.868*pow(10,24),
    5.683*pow(10,26)
    }; 

double Po[C][3] = { // posiciones iniciales
    {-1.218 * pow(10, 6), 6.361*pow(10,5), 2.326*pow(10,4)},
    {1.471 * pow(10, 8),-2.514*pow(10,7), 2.398*pow(10,4)},
    {-2.471 * pow(10, 8), -1.422*pow(10,7), 5.744*pow(10,6)},
    {6.427*pow(10,8),  -3.853*pow(10,8), -1.278*pow(10,7)},
    {7.909*pow(10,6),-6.716*pow(10,7),-6.354*pow(10,6)},
    {1.874*pow(10,7),-1.063*pow(10,8),-2.592*pow(10,6)},
    {9.780*pow(10,8),-1.119*pow(10,9),-1.948*pow(10,7)}
    };

double Vo[C][3] = { // velocidades iniciales
    {-8.131*pow(10,-3), -1.386*pow(10,-2), 3.036*pow(10,-4)},
    {4.598, 2.922*10, 1.65*pow(10,-5)},
    {2.354, -2.212*10, -5.211*pow(10,-1)},
    {6.560,  1.182*10, -1.958*pow(10,-1)},
    {3.849*pow(10,1),8.984,-2.796},
    {3.418*10,6.287,-1.886},
    {6.733,6.338,-3.788*pow(10,-1)}
    };


// CALCULA DISTANCIAS AL CUADRADO #####################################################################################

double dist2(double K[C][3], int a, int b){

    double dist2=0.;
    for (int j=0; j<3; j++){

        dist2+=pow((K[a][j]-K[b][j]),2); 

    }

    return dist2;

}

// CALCULA ACELERACION EN EL EJE INDICADO   ##################################################################
double a(int cuerpo, int eje,double L[C][3], double m[C]){ 

    double a=0.;

    for (int b=0; b<C; b++){

        if (cuerpo != b){ //no considerar el propio cuerpo

            a+=-m[b]*(L[cuerpo][eje]-L[b][eje])/pow(sqrt(dist2(L,cuerpo,b)),3);

        }
        
    }

    return a;
}

//CALCULA ENERGIAS############################################################################################
double Energia(double r[C][3], double v[C][3], double m[C], int cos){
    double E =0.;
    for (int j =0; j<3; j++){
        E += m[cos]*pow(v[cos][j],2);
    }
    
    for (int cos_2 =0; cos_2<C; cos_2++){
        if (cos_2 != cos){
            double d = pow(sqrt(dist2(r,cos,cos_2)), 3);
            E+= m[cos_2]*m[cos]/d;
        }
    }
    return E;
}
//####################################################################################################################3

int main(){
    FILE *output;
    output = fopen("trayectorias.txt", "w");

    double h = 2.0*PI*(1./(365.0*24*60));

    fprintf(output,"t ");

    for (int k=0; k<C; k++){ //header
    fprintf(output, "x_%d y_%d z_%d E_%d ", k,k,k,k );

    }
    
    fprintf(output,"\n");

    double P[C][3];
    double V[C][3];
    double M[C];

    for(int w=0; w<C; w++){ //condiciones iniciales normalizadas

        M[w]=MM[w]/MM[0];

        for(int z=0; z<3; z++){

            P[w][z]=Po[w][z]/L0;
            V[w][z]=Vo[w][z]/L0*T0;

        }

    }

    double t = 0.0;

    fprintf(output, "%lf ", t);

    for (int k = 0; k < C; k++){

        double E= Energia(P,V,M,k);
        fprintf(output, "%lf %lf %lf %lf ", P[k][0], P[k][1], P[k][2], E);

    }

    fprintf(output, "\n");

    double Kr[C][4][3];
    double Kv[C][4][3];
    double Pk[C][3];

    for (int n=0; n<N; n++){ //iterar por paso

        t+=h; //guarda el tiempo

        // CALCULAMOS K1, al mismo tiempo calculamos las "nuevas posiciones" para las siguientes k
       for (int k = 0; k < C; k++){ //iteramos por cuerpo
            for (int e=0; e<3; e++){ //iteramos por eje

                Kr[k][0][e]= V[k][e];
                Kv[k][0][e]= a(k,e,P,M);
        
                Pk[k][e]=P[k][e]+h/2*Kv[k][0][e];

            }
        } 

        // K2, idem antem
       for (int k = 0; k < C; k++){ //iteramos por cuerpo
            for ( int e=0; e<3; e++){ //iteramos por eje

                Kr[k][1][e]= V[k][e]+h/2*Kr[k][0][e];
                Kv[k][1][e]= a(k,e,Pk,M);
        
                Pk[k][e]=P[k][e]+h/2*Kv[k][1][e];

            }
        }

        // K3, idem antem
       for (int k = 0; k < C; k++){ //iteramos por cuerpo
            for (int e=0; e<3; e++){ //iteramos por eje

                Kr[k][2][e]= V[k][e]+h/2*Kr[k][1][e];
                Kv[k][2][e]= a(k,e,Pk,M);
        
                Pk[k][e]=P[k][e]+h*Kv[k][2][e];

            }
        }

        // K4, ultima k a calcular
       for (int k = 0; k < C; k++){ //iteramos por cuerpo
            for (int e=0; e<3; e++){ //iteramos por eje

                Kr[k][3][e]= V[k][e]+h*Kr[k][2][e];
                Kv[k][3][e]= a(k,e,Pk,M);

            }
        }

        //Una vez ya tenemos todas las k, calculamos y guardamos los resultados

        for (int k = 0; k < C; k++){ // iteramos por cuerpo
            for (int e=0; e<3; e++){  //iterar por eje para obtener nuevas posiciones y velocidades

                P[k][e] += h / 6 * (Kr[k][0][e] + 2*Kr[k][1][e] + 2 * Kr[k][2][e] + Kr[k][3][e]);
                V[k][e] += h/6*(Kv[k][0][e] + 2*Kv[k][1][e] + 2 * Kv[k][2][e] + Kv[k][3][e]);

            }
        }

        if(int(n%(24*60))==0){

            fprintf(output, "%lf ", t);

            for (int k=0; k<C; k++){
                    
                double E= Energia(P,V,M,k);
                fprintf(output, "%lf %lf %lf %lf ", P[k][0], P[k][1], P[k][2], E);

            }

            fprintf(output, "\n");

        }

    }

    return 0;

}


