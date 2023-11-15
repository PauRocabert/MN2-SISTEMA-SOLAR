#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>



const double PI = 3.14159265359;
//LONGITUD PRÒPIA: 1U.A
const double L0 = 149597870.7; //km
//TEMPS PROPI: 1 Any/2Pi
const double T0 = (365.0*24.0*3600.0)/(2*PI);//s 
int N_cossos = 4;
int N_dim =3;
char *names[] = {"Sol","Terra","Mart", "Jupiter"};



// ******************* FUNCIONS PRÈVIES ********************

double dvdt(double r[N_cossos][N_dim], double Distancies3[N_cossos][N_cossos], double m[N_cossos],int cos_central, int coord){
    double a_i = 0.0;
    double r1[N_dim];
    double r2[N_dim];
    for(int i =0; i<N_dim; i++) r1[i] =r[cos_central][i];
    for (int cos = 0; cos<N_cossos; cos++){
        if (cos != cos_central){  //No podem considerar la força que un cos fa sobre ell mateix
            for (int n=0; n<N_dim; n++) r2[n] = r[cos][n];
            a_i += m[cos]*(r2[coord] - r1[coord])/(Distancies3[cos_central][cos]);
            }
        }
    return a_i;
}
double Energia(double v[N_cossos][N_dim],double distancies[N_cossos][N_cossos], double m[N_cossos], int cos){
    double E =0.0;
    for (int j =0; j<N_dim; j++){
        E += m[cos]*pow(v[cos][j],2);
    }
    for (int cos_2 =0; cos_2<N_cossos; cos_2++){
        if (cos_2 != cos){
            E+= m[cos_2]*m[cos]/pow(distancies[cos][cos_2],1.0/3.0);
        }
    }
    return E;
}


int main(){
    clock_t t0 = clock();
//**DECLARACIÓ DE VARIABLES PRINCIPALS [POSICIONS, VELOCITATS I MASSES NORMALITZADES]**
    double r[N_cossos][N_dim];
    double v[N_cossos][N_dim];
    double E0[N_cossos];
    double m[N_cossos]; 

//***************************DADES + CONDICIONS INICIALS***************************
    double M[] = {1.99*pow(10.,30.), 5.97*pow(10.,24.), 6.42*pow(10.,23.), 1.9*pow(10.,27.)};
    double R0[3][4] = {
        {-1.218*pow(10,6), 1.471*pow(10,8), -2.471*pow(10,8), 6.427*pow(10,8)},
        {6.361*pow(10,5), -2.514*pow(10,7), -1.422*pow(10,7), -3.853*pow(10,8)},
        {2.326*pow(10.0,4.0), 2.398*pow(10.0, 4.0),  5.744*pow(10.0,6.0), -1.2778*pow(10.0,7.0)}};
    double V0[3][4] ={
        {-8.131*pow(10,-3), 4.598, 2.354, 6.560},
        {-1.386*pow(10,-2), 2.922*10, -2.212*10, 1.182*10},
        {3.036*pow(10,-4), 1.65*pow(10,-5), -5.211*pow(10,-1), -1.958*pow(10,-1)}};

//Normalització de les masses
    for (int i=0; i<N_cossos; i++) m[i] =M[i]/M[0]; 
//Normalització de les posicions i les velocitats
    for (int j=0; j<N_dim; j++){
        for (int cos=0; cos<N_cossos; cos++){
        r[cos][j] = R0[j][cos]/L0;
        v[cos][j] = V0[j][cos]/L0*T0;
    }}

//*********DISCRETITZACIÓ********** 

//    double dt =  2.0*PI*(1./365.0); //dt = 1dia (expressat en funció de T0)
    double dt =  2.0*PI*(1./(365.0*24.0*60.0)); //dt = 1min (expressat en funció de T0)
//    double dt = 2.0*PI*(1./365.0)*(1/(24.0*3600.0)); //dt = 1s (expressat en funció de T0)
    int N = 782*24*60; // Minuts entre  entre 13/9/2021 i 3/11/2003



    double distanceM[N_cossos][N_cossos];
    for (int cos1 =0; cos1< N_cossos; cos1++){
        for (int cos2=0; cos2<N_cossos; cos2++) distanceM[cos1][cos2] = .0;} 

    void DistanceMatrix(double r[N_cossos][N_dim]){
        for (int cos1 = 0; cos1 <N_cossos; cos1++){
            for (int cos2=cos1+1; cos2<N_cossos; cos2++){
                double d = 0.;
                for (int i=0; i<N_dim; i++){
                    d+= pow((r[cos1][i] - r[cos2][i]),2); 
                }
                distanceM[cos1][cos2] = pow(d,1.5);
                distanceM[cos2][cos1] = distanceM[cos1][cos2];
            }}
    };

//***********DOCUMENT - CAPÇALERA**************
    FILE *pF = fopen("sistema_solar_posicions_postopt.txt","w");
    fprintf(pF, "t ");
    for (int cos =0; cos<N_cossos; cos++){
        fprintf(pF,"x_%s y_%s z_%s E_%s difE_%s ", names[cos], names[cos], names[cos], names[cos], names[cos]);
    }
    fprintf(pF, "\n%lf ", 0.0);
    DistanceMatrix(r);
    for (int cos =0; cos<N_cossos; cos++){
        for (int j=0; j<N_dim; j++){
            fprintf(pF, "%lf ", r[cos][j]);
        }
        E0[cos] = Energia(v, distanceM, m, cos);
        fprintf(pF, "%lf %lf ", E0[cos], E0[cos]-E0[cos]);
    }
//**************RUNGE KUTTA************
    double K1[N_cossos][N_dim][2];
    double K2[N_cossos][N_dim][2];
    double K3[N_cossos][N_dim][2];
    double K4[N_cossos][N_dim][2];
    double rK[N_cossos][N_dim];

    for (int i=1; i<N+1; i++){
        double t = i*dt; 
        DistanceMatrix(r);
        for (int cos=0; cos<N_cossos; cos++){
            for (int j=0; j<N_dim; j++){
                K1[cos][j][0] = v[cos][j];
                K1[cos][j][1] = dvdt(r, distanceM, m, cos, j);
                rK[cos][j] = r[cos][j] +  0.5*dt*K1[cos][j][1];
                }}
        DistanceMatrix(rK);
        for (int cos=0; cos<N_cossos; cos++){
            for (int j=0; j<N_dim; j++){
                K2[cos][j][0] = v[cos][j]+K1[cos][j][1]*dt/2;
                K2[cos][j][1] = dvdt(rK, distanceM, m, cos, j);
                rK[cos][j] = r[cos][j] +  0.5*dt*K2[cos][j][1];
                }}
        DistanceMatrix(rK);
        for (int cos=0; cos<N_cossos; cos++){
            for (int j=0; j<N_dim; j++){
                K3[cos][j][0] = v[cos][j]+K2[cos][j][1]*dt/2;
                K3[cos][j][1] = dvdt(rK,distanceM, m, cos, j);
                rK[cos][j] = r[cos][j] +  0.5*dt*K3[cos][j][1];
                }}

        DistanceMatrix(rK);
        for (int cos=0; cos<N_cossos; cos++){
            for (int j=0; j<N_dim; j++){
                K4[cos][j][0] = v[cos][j]+K3[cos][j][1]*dt;
                K4[cos][j][1] = dvdt(rK, distanceM, m, cos, j);}} 

        for (int cos=0; cos<N_cossos; cos++){
            for (int j=0; j<N_dim; j++){
                r[cos][j] += (dt/6)*(K1[cos][j][0]+2*K2[cos][j][0]+2*K3[cos][j][0]+K4[cos][j][0]);
                v[cos][j] += (dt/6)*(K1[cos][j][1]+2*K2[cos][j][1]+2*K3[cos][j][1]+K4[cos][j][1]);
            }}
        
        if (i%1440 ==0){
            fprintf(pF,"\n%lf ",t);
            for (int cos=0; cos<N_cossos; cos++){
                for (int j=0; j<N_dim; j++){
                fprintf(pF,"%lf ", r[cos][j]);
                }
            fprintf(pF,"%lf %lf ", Energia(v, distanceM, m, cos), Energia(v, distanceM, m, cos)-E0[cos]);
            }
        }
    }
    t0 = clock()-t0;
    double temps_exc = ((double)t0)/CLOCKS_PER_SEC;
    printf("Temps execucio t=%lf", temps_exc);
    return 0;
}

