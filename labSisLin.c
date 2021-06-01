#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

int main(){

    SistLinear_t *SL;

    while(1){
        double *t;
        SL = lerSistLinear();
        real_t *X = alocaVetor(SL->n);

        printf("SISTEMA:\n");
        prnSistLinear(SL);

        if(eliminacaoGauss(SL, X, t)){
            perror("erro eliminacao gaus");
        }

        printf("GAUS:\n");
        prnSistLinear(SL);

        printf("X:\n");
        prnVetor(X, SL->n);

        liberaSistLinear(SL);
        liberaVetor(X);
    }
}