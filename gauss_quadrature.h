#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

#include <stdlib.h>
#include <gsl/gsl_integration.h>

typedef struct {
    double *points;
    double *weights;
} GaussQuadrature;

// Funktion zur Berechnung von Gauss-Quadraturpunkten und -gewichten
GaussQuadrature calculate_gauss_legendre_points_weights(int n, double a, double b) {
    GaussQuadrature result;
    result.points = (double *)malloc(n * sizeof(double));
    result.weights = (double *)malloc(n * sizeof(double));
    
    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(n);
    
    for (int i = 0; i < n; i++) {
        double xi, wi;
        gsl_integration_glfixed_point(a, b, i, &xi, &wi, t);
        result.points[i] = xi;
        result.weights[i] = wi;
    }

    gsl_integration_glfixed_table_free(t);
    return result;
}

// Funktion zum Freigeben des Speichers
void free_gauss_quadrature(GaussQuadrature *gq) {
    free(gq->points);
    free(gq->weights);
}

#endif // GAUSS_QUADRATURE_H
