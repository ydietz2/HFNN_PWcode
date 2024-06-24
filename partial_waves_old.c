#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <time.h>
#include "gauss_quadrature.h"

int Nrows = 62;

void read(const char *filepath, double matrix2[Nrows * Nrows][3] ){
    FILE *file = fopen(filepath, "r"); 

    if (file == NULL) {
        perror("Fehler beim Öffnen der Datei");
        return ;
    }


    int rows2 = Nrows * Nrows; 
    int cols2 = 3; 

    // Skipping Rows (Mesh)

    for (int skipper = 0; skipper < Nrows; skipper++){
        fscanf(file, "%*[^\n]\n");
    }


    // Reading File as matrix
    for (int i = 0; i < rows2; i++) {
        for (int j = 0; j < cols2; j++) {
            fscanf(file, "%lf", &matrix2[i][j]);
        }
    }

    fclose(file);
}

// as the name says only read diag entries 
void readdiag(const char *filepath, double x[Nrows], double y[Nrows] ){
    FILE *file = fopen(filepath, "r"); 
    double matrix2[Nrows * Nrows][3];

    if (file == NULL) {
        perror("Fehler beim Öffnen der Datei");
        return ;
    }


    int rows2 = Nrows * Nrows; 
    int cols2 = 3; 

    // Skipping Rows (Mesh)

    for (int skipper = 0; skipper < Nrows; skipper++){
        fscanf(file, "%*[^\n]\n");
    }


    // Reading File as matrix
    for (int i = 0; i < rows2; i++) {
        for (int j = 0; j < cols2; j++) {
            fscanf(file, "%lf", &matrix2[i][j]);
        }
    }

    fclose(file);

    for (int i = 0; i*Nrows < rows2; i++){
        
        x[i] = matrix2[i*Nrows+i][0];
        y[i] = matrix2[i*Nrows+i][2];

    }
    x[0] = 0.0;
}


// x and y from VNN, N = Rows, len = lenght of xeval, xeval = points to evaluate interpolation, yfinal = results
void interpolate (double x[], double y[], int N, int len ,double xeval[], double yfinal[]){

    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);

    gsl_spline_init (spline, x, y, N);

    double xi;

    for (int i = 0; i < len; i++){
        yfinal[i] = gsl_spline_eval (spline, xeval[i], acc);
    }

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

}


//Meshfunction
void mesh(float a, float b, int d, double xi[d], double wi[d]){

    
    GaussQuadrature gq = calculate_gauss_legendre_points_weights(d, a, b);

    for (int i = 0; i < d; i++){
        wi[i] = gq.weights[i] ;
        xi[i] =  gq.points[i] ;
    }
    free_gauss_quadrature(&gq);    
}


int main() {
    clock_t begin = clock();
    const double PI = 3.1415926;
    double unit = 197.326 * 197.326 / ((938.272 + 939.565) / 2.0); 
    // File specific
    int cols = 3;
    int rows = Nrows * Nrows;
    
    // gridsize 
    int grid = 20;
    int wgrid = 20;
    int counter = 0;

    // p def
    double pi[grid];
    double pwi[grid];
    int lpi = sizeof(pi) / sizeof(pi[0]); // lenght of p array 

    // k def
    double ki[grid];
    double kwi[grid];
    int lxi = sizeof(ki) / sizeof(ki[0]);

    //theta def 
    double ti[wgrid];
    double twi[wgrid];
    int lti = sizeof(ti) / sizeof(ti[0]);

    //quantum numbers
    int S = 0;
    int L = 0;
    int Lprime = 0;
    //int J =0;
    int T = 1;
    double VNN[Nrows];
    char str [10000];
    

    //interpolation arrys
    double x[Nrows];
    double y[Nrows];
    double fki[grid]; 

    // Save File
    FILE *savefile;
    savefile = fopen("VLO_HF_energy_pw.txt","w");

    int Jmax = 8;

    float kf;
    //float kfar[] = {0.66651051, 0.83975062, 0.96127449, 1.05801948, 1.13971693,
    //  1.21112997, 1.27498873, 1.33302101, 1.38639772, 1.43595336,
    //   1.4823061 , 1.52592814, 1.56718927, 1.60638514, 1.64375626,
    //   1.67950124, 1.7137862 , 1.74675168, 1.77851774, 1.80918786};
    float kfar[] = {1.05801948, 1.2111299, 1.33302101};
    for (int runner= 0; runner < 3; runner++){
    kf = kfar[runner];
    double rho = 2.0 * kf*kf*kf / (3.0 * PI * PI);
    float Areagauss = 0.0;
    mesh(0.0, 2.0*kf, grid, pi, pwi); //p mesh
    for (int j=0; j < lpi; j++){ // p loop
        float kmax = sqrt(kf*kf - pi[j]* pi[j] / 4);
        mesh(0, kmax, grid, ki, kwi); //k mesh
        float prog = (double)j / 10.0;
        printf("Progress: %f \n", prog);
        for (int i=0; i < lxi; i++){ // k loop
            float ctheta = (kf*kf - pi[j]*pi[j] / 4 - ki[i]* ki[i] ) / (ki[i]* pi[j]);  
            //double xmin = fmax(-1.0, -ctheta);
            //double xmax = fmin(1.0, ctheta); 
            double xmin = fmax(-1.0, fmin(1.0,-ctheta));
            double xmax = fmin(1.0, fmax(-1.0,ctheta)); 
            double max_min_diff = xmax - xmin;
            if (xmax < xmin){
                continue;
            }
            mesh(xmin, xmax, wgrid,ti, twi );
            //for (int k = 0; k < lti; k++ ){
                for (int qJ = 0; qJ <= Jmax; qJ++){
                    for (int qS = 0; qS < 2; qS++){
                        for (int qT = 0; qT <2; qT++){
                              if (qS == 0){
                                int qL = qJ;
                                Lprime = qL;
                                if ((qT + qL )%2 == 0){
                                    continue;
                                }

                                double factor = 1.0;
                                //if (qT==0) {
                                //    factor = -3.0;
                                //}
                                sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_np.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                readdiag(str, x,y);
                                interpolate(x,y,Nrows,grid,ki,fki);
                                //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                Areagauss += (1 / (4* PI * PI *PI)) * factor * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 
    
                                if (qT == 1){
                                    sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_nn.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                    readdiag(str, x,y);
                                    interpolate(x,y,Nrows,grid,ki,fki);
                                    //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                    Areagauss += (1 / (4* PI * PI *PI)) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 

                                    sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_pp.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                    readdiag(str, x,y);
                                    interpolate(x,y,Nrows,grid,ki,fki);
                                    //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                    Areagauss += (1 / (4* PI * PI *PI)) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 
                                }
                              }
                              else{
                              
                                for (int steps = -1 ; steps < 2; steps+=1){
                                int qL = qJ + steps;
                                Lprime = qL;
                                if (qL < 0){
                                    continue;
                                }
                                if (qL == 0 && qJ == 0){
                                    continue;
                                    }

                                if ((qT + qL + qS )%2 == 0){
                                    continue;
                                }

                                double factor = 1.0;
                                //if (qT==0) {
                                //    factor = -3.0;
                                //}

                                sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_np.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                readdiag(str, x,y);
                                interpolate(x,y,Nrows,grid,ki,fki);
                                //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                Areagauss += (1 / (4* PI * PI *PI)) * factor * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 
    
                                if (qT == 1){
                                    sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_nn.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                    readdiag(str, x,y);
                                    interpolate(x,y,Nrows,grid,ki,fki);
                                    //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                    Areagauss += (1 / (4* PI * PI *PI)) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 

                                    sprintf(str,"magic_ME/VNN_N3LO_EM500_SLLJT_%d%d%d%d%d_lambda_1.80_Np_%d_pp.dat", qS, qL, Lprime, qJ, qT, Nrows);
                                    readdiag(str, x,y);
                                    interpolate(x,y,Nrows,grid,ki,fki);
                                    //Areagauss += 1 / (4* PI * PI *PI) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * twi[k] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ;
                                    Areagauss += (1 / (4* PI * PI *PI)) * ki[i]* ki[i] * pi[j]* pi[j] * kwi[i] * pwi[j] * max_min_diff * (2.0 * qJ +1.0) * fki[i] * unit * 2.0 ; 
                                }
                              }



                            }
                        }   
                    }

                 }
           // }
        }
    }
    printf("%lf %lf \n", kf, Areagauss / rho);
    printf("density = %lf \n", rho);
    fprintf(savefile,"%lf %lf \n", kf, Areagauss/rho );
    }
    fclose(savefile);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time spent: %lf \n", time_spent ); // 1365s 
    }
