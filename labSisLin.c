#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

int main(){

    SistLinear_t *SL, *SL2;
    real_t *R;
    int i;
    for(i = 1; 1 ; i++){
        double *t;
        SL = lerSistLinear();
        real_t L2, *X = alocaVetor(SL->n);

        printf("***** Sistema %d --> n = %d, erro: %f\n", i, SL->n, SL->erro);

        //copio a sistema linear para manter os meus coeficientes originais
        SL2 = copiaMatriz(SL);
        if(eliminacaoGauss(SL2, X, t)){
            perror("erro eliminacao gaus");
        }

        printf("===> Eliminação Gauss: %f ms\n--> X: ", 0.0);
        prnVetor(X, SL2->n);
        L2 = normaL2Residuo(SL2, X, R);
        //if(L2 < 5)
        printf("--> Norma L2 do residuo: %f\n", L2);
        pulaLinha(1);

        gaussJacobi(SL, X, t);
        printf("===> Jacobi: %f ms\n--> X: ", 0.0);
        prnVetor(X, SL2->n);
        L2 = normaL2Residuo(SL2, X, R);
        printf("--> Norma L2 do residuo: %f\n", L2);
        pulaLinha(1);

        gaussSeidel(SL, X, t);
        printf("===> Gauss-Seidel: %f ms\n--> X: ", 0.0);
        prnVetor(X, SL2->n);
        L2 = normaL2Residuo(SL2, X, R);
        printf("--> Norma L2 do residuo: %f\n", L2);
        pulaLinha(1);

        liberaSistLinear(SL2);
        liberaSistLinear(SL);
        liberaVetor(X);
    }
}