#include <stdio.h>
#include <time.h>
#include <math.h>

//DEFINICIÓ DEL SISTEMA 

#define cossos 7 //Nombre de cossos
#define dim 3 //Dimensions

//CONSTANTS

const double PI = 3.14159265359;
const double L0 = 149597870.7; //km //Longitud pròpia: 1U.A
const double T0 = (365.2425*24.*3600.)/(2*PI);//s  //Temps propi: 1 Any/2Pi

double dt = 2*PI/(365.2425*24*60); //min 
#define N 782*24*60//Nombre d'iteracions (nombre de min en 782 dies)

/*
double dt = 3600./T0 //hora
#define N 18768//Nombre d'iteracions (nombre de h en 782 dies)

double dt = 3600.*24./T0; //dia
#define N 782 //Nombre d'iteracions (nombre de dies en 782 dies)
*/

char *names[] = {"Sol","Mercuri","Venus","Terra","Mart", "Jupiter","Saturn"};

//DADES INICIALS sense normalitzar
double M[cossos]= { //masas
    1.99*pow(10.,30.),
    3.302*pow(10.,23.),
    4.868*pow(10.,24.),
    5.97*pow(10.,24.),
    6.42*pow(10.,23.),
    1.90*pow(10.,27.),
    5.683*pow(10.,26)
    }; 

double Po[cossos][dim] = { // posiciones iniciales
    {-1.218 * pow(10., 6), 6.361*pow(10.,5), 2.326*pow(10.,4)},
    {7.909*pow(10.,6),-6.716*pow(10.,7),-6.354*pow(10.,6)},
    {1.874*pow(10.,7),-1.063*pow(10.,8),-2.592*pow(10.,6)},
    {1.471 * pow(10., 8),-2.514*pow(10.,7), 2.398*pow(10.,4)},
    {-2.471 * pow(10., 8), -1.422*pow(10.,7), 5.744*pow(10.,6)},
    {6.427*pow(10.,8),  -3.853*pow(10.,8), -1.278*pow(10.,7)},
    {9.780*pow(10.,8),-1.119*pow(10.,9),-1.948*pow(10.,7)}
    };

double Vo[cossos][dim] = { // velocidades iniciales
    {-8.131*pow(10.,-3), -1.386*pow(10.,-2), 3.036*pow(10.,-4)},
    {3.849*pow(10.,1),8.984,-2.796},
    {3.418*10,6.287,-1.886},
    {4.598, 2.922*10, 1.65*pow(10.,-5)},
    {2.354, -2.212*10, -5.211*pow(10.,-1)},
    {6.560,  1.182*10, -1.958*pow(10.,-1)},
    {6.733,6.338,-3.788*pow(10.,-1)}
    };

double a(double r[cossos][dim], double Distancies3[cossos][cossos], double m[cossos],int k1, int eix){
    double a_i = 0.0;
    for (int k2 = 0; k2<cossos; k2++){ 
        if (k2 != k1){  //No podem considerar la força que un cos fa sobre ell mateix
            a_i += m[k2]*(r[k2][eix] - r[k1][eix])/(Distancies3[k1][k2]); //distàncies al cub
            }
        }
    return a_i;}


double Energia(double v[cossos][dim],double distancies3[cossos][cossos], double m[cossos], int k1){
    double E =0.0;
    for (int e =0; e<dim; e++){
        E += 0.5*m[k1]*pow(v[k1][e],2);
    }
    for (int k2 =0; k2<cossos; k2++){
        if (k2 != k1){
            E+= m[k2]*m[k1]/pow(distancies3[k1][k2],1.0/3.0);
        }
    }
    return E;
}



int main(){
    clock_t t0 = clock();
    FILE *output;
    output = fopen("trajectories_difE.txt", "w");

    fprintf(output,"t");
    for (int k=0; k<cossos; k++){ //headers
    fprintf(output, " x_%s y_%s z_%s E_%s", names[k],names[k],names[k],names[k]);
    }
    fprintf(output,"\n");


    double m[cossos];
    double P[cossos][dim];
    double V[cossos][dim];
    
    //Normalització de les masses
    for (int k=0; k<cossos; k++) m[k] =M[k]/M[0]; 
    
    //Normalització de les posicions i les velocitats
    for (int e=0; e<dim; e++){
        for (int k=0; k<cossos; k++){
            P[k][e] = Po[k][e]/L0;
            V[k][e] = Vo[k][e]/L0*T0;
        }
    }
    
    //MATRIU DE DISTÀNCIES
    double dist[cossos][cossos];
    for (int k1 =0; k1< cossos; k1++){
        for (int k2=0; k2<cossos; k2++) dist[k1][k2] = .0;} 

    void MDIST3(double r[cossos][dim]){
        for (int k1 = 0; k1 <cossos; k1++){
            for (int k2=k1+1; k2<cossos; k2++){
                double d = 0.;
                for (int i=0; i<dim; i++){
                    d+= pow((r[k1][i] - r[k2][i]),2); 
                }
                dist[k1][k2] = pow(d,1.5);
                dist[k2][k1] = dist[k1][k2];
            }
        }
    };
    double t = 0.;

    fprintf(output, "%lf", t);

    MDIST3(P);
    for (int k = 0; k < cossos; k++){
        double E0= Energia(V,dist,m,k);
        fprintf(output, " %lf %lf %lf %lf", P[k][0], P[k][1], P[k][2], E0);
    }

    fprintf(output, "\n");


    //RK4
    double Kp[cossos][4][dim]; //Kp[cos][subíndex K][dimensió]
    double Kv[cossos][4][dim]; //Kv[cos][subíndex K][dimensió]
    double Pk[3][cossos][dim];

    for (int n=1; n<N; n++){ //iterar per pas
        t+=dt; //guarda el temps

       MDIST3(P);
        // CALCULEM K1, al mateix temps calculem les "noves posicions" per a les següents
       for (int k = 0; k < cossos; k++){ //iterem per cos
            for (int e=0; e<dim; e++){ //iterem per eix
                Kp[k][0][e]= V[k][e];
                Kv[k][0][e]= a(P, dist, m, k, e);
                Pk[0][k][e]=P[k][e]+dt/2*Kv[k][0][e];
            }
        } 
        // K2, idem antem
       MDIST3(Pk[0]);
       for (int k = 0; k < cossos; k++){ //iterem per cos
            for ( int e=0; e<dim; e++){ //iterem per eix
                Kp[k][1][e]= V[k][e]+dt/2*Kp[k][0][e];
                Kv[k][1][e]= a(Pk[0], dist, m,k,e);
                Pk[1][k][e]=P[k][e]+dt/2*Kv[k][1][e];
            }
        }
       MDIST3(Pk[1]);
        // K3, idem antem
       for (int k = 0; k < cossos; k++){ //iterem per cos
            for (int e=0; e<dim; e++){ //iterem per eix
                Kp[k][2][e]= V[k][e]+dt/2*Kp[k][1][e];
                Kv[k][2][e]= a(Pk[1], dist, m,k,e);
                Pk[2][k][e]=P[k][e]+dt*Kv[k][2][e];
            }
        }
        // K4, ultima k a calcular

       MDIST3(Pk[2]);
       for (int k = 0; k < cossos; k++){ //iterem per cos
            for (int e=0; e<3; e++){ //iterem per eix
                Kp[k][3][e]= V[k][e]+dt*Kp[k][2][e];
                Kv[k][3][e]= a(Pk[2], dist, m,k,e);
            }
        }

        //Una volta ja tenim totes les k, calculem i guardem els resultats
    
        for (int k = 0; k < cossos; k++){ // iteramos por cuerpo
            for (int e=0; e<3; e++){  //iterar por eje para obtener nuevas posiciones y velocidades
                P[k][e] += dt/6 * (Kp[k][0][e] + 2*Kp[k][1][e] + 2 * Kp[k][2][e] + Kp[k][3][e]);
                V[k][e] += dt/6*(Kv[k][0][e] + 2*Kv[k][1][e] + 2 * Kv[k][2][e] + Kv[k][3][e]);
            }
        }

        MDIST3(P);
        if((n%(24*60))==0){
            fprintf(output, "%lf", t);
            for (int k=0; k<cossos; k++){           
                double E= Energia(V,dist,m,k);
                fprintf(output, " %lf %lf %lf %lf", P[k][0], P[k][1], P[k][2], E);
            }
            fprintf(output, "\n");
        }
    }
    t0 = clock()-t0;
    double temps_exc = ((double)t0)/CLOCKS_PER_SEC;
    printf("Temps execucio t=%lf", temps_exc);
    return 0;
     
    
}