#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

int main(){

    SistLinear_t *SL, *SL2;
    real_t *R;
    int i, k;
    double t;
    for(i = 1; 1 ; i++){
        SL = lerSistLinear();
        real_t L2, *X = alocaVetor(SL->n);

        printf("***** Sistema %d --> n = %d, erro: %f\n", i, SL->n, SL->erro);

        //copio a sistema linear para manter os meus coeficientes originais
        SL2 = copiaMatriz(SL);

        if(eliminacaoGauss(SL2, X, &t)){
            perror("erro eliminacao gaus");
        }
        else{
            printf("===> Eliminação Gauss: %f ms\n--> X: ", t);
            prnVetor(X, SL2->n);
            R = residuo(SL, X);
            L2 = normaL2Residuo(SL, X, R);
            printf("--> Norma L2 do residuo: %1.8e\n", L2);
            pulaLinha(1);
        }

        if(L2 > 5.0){
            chamaRefinamento(SL, X, &t);
        }

        k = gaussJacobi(SL, X, &t);
        if(k != -1){
            printf("===> Jacobi: %f ms --> %d iteracoes\n--> X: ", t, k);
            prnVetor(X, SL->n);
            R = residuo(SL, X);
            L2 = normaL2Residuo(SL, X, R);
            printf("--> Norma L2 do residuo: %1.8e\n", L2);
            pulaLinha(1);
        }

        if(L2 > 5.0){
            chamaRefinamento(SL, X, &t);
        }

        k = gaussSeidel(SL, X, &t);
        if(k != -1){
            printf("===> Gauss-Seidel: %f ms --> %d iteracoes\n--> X: ", t, k);
            prnVetor(X, SL->n);
            R = residuo(SL, X);
            L2 = normaL2Residuo(SL, X, R);
            printf("--> Norma L2 do residuo: %1.8e\n", L2);
            pulaLinha(1);
        }

        if(L2 > 5.0){
            chamaRefinamento(SL, X, &t);
        }

        liberaSistLinear(SL2);
        liberaSistLinear(SL);
        liberaVetor(X);
    }
}