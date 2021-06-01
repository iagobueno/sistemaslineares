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

        printf("SISTEMA:\n");
        prnSistLinear(SL);

        SL2 = copiaMatriz(SL);

        if(eliminacaoGauss(SL2, X, t)){
            perror("erro eliminacao gaus");
        }

        printf("MATRIZ 2:\n");
        prnSistLinear(SL2);

        printf("X:\n");
        prnVetor(X, SL->n);

        liberaSistLinear(SL);
        liberaVetor(X);
    }
}