#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

int main(){

    SistLinear_t *SL, *SL2;

    while(1){
        double *t;
        SL = lerSistLinear();
        real_t *X = alocaVetor(SL->n);

        // printf("SISTEMA:");
        // pulaLinha(2);
        // prnSistLinear(SL);
        // pulaLinha(1);

        SL2 = copiaMatriz(SL);

        if(eliminacaoGauss(SL2, X, t)){
            perror("erro eliminacao gaus");
        }

        // printf("MATRIZ GAUSS:");
        // pulaLinha(2);
        // prnSistLinear(SL2);
        // pulaLinha(1);

        printf("GAUSS JORDAN:");
        pulaLinha(1);
        prnVetor(X, SL2->n);
        pulaLinha(2);

        gaussJacobi(SL, X, t);
        // printf("GAUSS JACOBI:");
        // pulaLinha(1);
        // prnVetor(X, SL2->n);


        gaussSeidel(SL, X, t);
        printf("GAUSS SEIDEL:");
        pulaLinha(1);
        prnVetor(X, SL2->n);

        liberaSistLinear(SL2);
        liberaSistLinear(SL);
        liberaVetor(X);
    }
}