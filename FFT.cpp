#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

const double pi = 3.14159265359;
const double L  = 2.*pi;
const int    N  = 8;
const double dx = L/(N );

void poisson(double grid_mass[][N][N], double grid_acc[][N][N][3]){

    int i, j, k;
    double *X, *Y, *Z ;
    X = (double*) malloc(N*sizeof(double));
    Y = (double*) malloc(N*sizeof(double));
    Z = (double*) malloc(N*sizeof(double));
    for(i = 0; i < N; i++){
        X[i] = -pi + i*dx ;
        for(j = 0; j < N; j++){
            Y[j] = -pi + j*dx ;
            for (k = 0; k < N; k++){
                Z[k]= -pi + k*dx ;
            }
        }
    }
    
    fftw_complex  *out1, *in2, *out2, *in1;
    fftw_plan     p1, p2;

    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N*N*N) );
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N*N*N) );
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N*N*N) );
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N*N*N) );

    p1 = fftw_plan_dft_3d(N, N, N, in1, out1, FFTW_FORWARD,FFTW_MEASURE );
    p2 = fftw_plan_dft_3d(N, N, N, in2, out2, FFTW_BACKWARD,FFTW_MEASURE);

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            for (k = 0; k < N; k++){
                in1[i*N*N + j*N + k][0] = grid_mass[i][j][k]; // row major ordering
                in1[i*N*N + j*N + k][1] = 0 ;
            }
        }
    }

    fftw_execute(p1); // FFT forward

    for ( i = 0; i < N; i++){   // f = g / ( kx^2 + ky^2 + kz^2 )
        for( j = 0; j < N; j++){
            for ( k = 0; k < N; k++){
                double ksquared=0;
                in2[i*N*N + j*N + k][0]=0;
                in2[i*N*N + j*N + k][1]=0;
                if(2*i<N){
                    ksquared=((double)i*i);
                }else{
                    ksquared=((double)(N-i)*(N-i));
                }
                if(2*j<N){
                    ksquared+=((double)j*j);
                }else{
                    ksquared+=((double)(N-j)*(N-j));
                }
                if (2*k<N){
                    ksquared+=((double)k*k);
                }else{
                    ksquared+=((double)(N-k)*(N-k));
                }
                if(ksquared!=0){
                    in2[i*N*N + j*N + k][0] = out1[i*N*N + j*N + k][0]/ksquared;
                    in2[i*N*N + j*N + k][1] = out1[i*N*N + j*N + k][1]/ksquared;
                }else{
                    in2[i*N*N + j*N + k][0] = 0;
                    in2[i*N*N + j*N + k][1] = 0;
                }
            }
        }
    }

    fftw_execute(p2); //FFT backward

    // checking the results computed

    double erl1 = 0.;
    for ( i = 0; i < N; i++) {
        for( j = 0; j < N; j++){
            for ( k = 0; k < N; k++){
                erl1 += fabs( in1[i*N*N + j*N + k][0] -  3.*out2[i*N*N + j*N + k][0]/(N*N*N))*dx*dx;
                std::cout<< i <<" "<< j<<" "<< k<<" "<< sin(X[i])*sin(Y[j])*cos(Z[k])<<" "<<  out2[i*N*N+j*N+k][0]/(N*N*N) <<" "<< std::endl; // > output
            }
        }
    }
    printf("\n");
    std::cout<< "erl1 = "<< erl1 << std::endl ;  // L1 error
    printf("-----------------------------------------------------------------------------------\n");
    
    double phi[N][N][N];
    for ( i = 0; i < N; i++){
        for ( j = 0; j < N; j++){
            for ( k = 0; k < N; k++){
                phi[i][j][k]=out2[i*N*N+j*N+k][0]/(N*N*N);
            }
        }
    }

    for ( i = 0; i < N; i++){
        for ( j = 0; j < N; j++){
            for ( k = 0; k < N; k++){
                if (i==0){
                    grid_acc[i][j][k][0]=-(phi[1][j][k]-phi[N-1][j][k])/(2.*dx);
                }else if (i==(N-1)){
                    grid_acc[i][j][k][0]=-(phi[0][j][k]-phi[N-2][j][k])/(2.*dx);
                }else{
                    grid_acc[i][j][k][0]=-(phi[i+1][j][k]-phi[i-1][j][k])/(2.*dx);
                }
                
                if (j==0){
                    grid_acc[i][j][k][1]=-(phi[i][1][k]-phi[i][N-1][k])/(2.*dx);
                }else if (j==(N-1)){
                    grid_acc[i][j][k][1]=-(phi[i][0][k]-phi[i][N-2][k])/(2.*dx);
                }else{
                    grid_acc[i][j][k][1]=-(phi[i][j+1][k]-phi[i][j-1][k])/(2.*dx);
                }
                
                if (k==0){
                    grid_acc[i][j][k][2]=-(phi[i][j][1]-phi[i][j][N-1])/(2.*dx);
                }else if (i==(N-1)){
                    grid_acc[i][j][k][2]=-(phi[i][j][0]-phi[i][j][N-2])/(2.*dx);
                }else{
                    grid_acc[i][j][k][2]=-(phi[i][j][k+1]-phi[i][j][k-1])/(2.*dx);
                }
            }
        }
    }
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(out1);
    fftw_free(out2);
    fftw_free(in1);
    fftw_free(in2);

}

int main(){
    double grid_mass[N][N][N], grid_acc[N][N][N][3];
    
    int i, j, k;
    double *X, *Y, *Z ;
    X = (double*) malloc(N*sizeof(double));
    Y = (double*) malloc(N*sizeof(double));
    Z = (double*) malloc(N*sizeof(double));
    for(i = 0; i < N; i++){
        X[i] = -pi + i*dx ;
        for(j = 0; j < N; j++){
            Y[j] = -pi + j*dx ;
            for (k = 0; k < N; k++){
                Z[k]= -pi + k*dx ;
                grid_mass[i][j][k] = 3.*sin(X[i])*sin(Y[j])*cos(Z[k]);
            }
        }
    }
    
    poisson(grid_mass,grid_acc);
    double err(0);
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                std::cout<<i<<j<<k
                <<" " <<grid_acc[i][j][k][0]<<" "<<grid_acc[i][j][k][1]<<" "<<grid_acc[i][j][k][2]
                <<" "<< -(sin(X[i+1])*sin(Y[j])*cos(Z[k])-sin(X[i-1])*sin(Y[j])*cos(Z[k]))/(2.*dx)<<" "<< -(sin(X[i])*sin(Y[j+1])*cos(Z[k])-sin(X[i])*sin(Y[j-1])*cos(Z[k]))/(2.*dx) <<" "<< -(sin(X[i])*sin(Y[j])*cos(Z[k+1])-sin(X[i])*sin(Y[j])*cos(Z[k-1]))/(2.*dx)<<std::endl;
            }
        }
    }
    return 0;
}
