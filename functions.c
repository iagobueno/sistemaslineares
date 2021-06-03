#include "SistemasLineares.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//testa se a alocacao dinamica deu certo
void testaMalloc(void *t){
    if(t == NULL){
        perror("alocacao de memoria invalida");
        exit(2);
    }
}

//aloca os coeficientes de uma matriz com a
//estratégia de vetor de ponteiros de linhas contíguas
real_t **alocaCoef(int n){
    real_t **c;    
    c = malloc( n * sizeof(real_t*));
    testaMalloc(c);
    c[0] = malloc( n * n * sizeof(real_t));
    testaMalloc(c[0]);

    int i, j;
    for(i = 0;i < n;i++)
        c[i] = c[0] + n * i;

    return c;
}

//aloca um vetor
real_t *alocaVetor(int n){
    real_t *vetor;
    
    vetor = malloc( n * sizeof(real_t));
    testaMalloc(vetor);

    return vetor;
}

//encontra o maior termo de uma coluna
int encontraMaior(SistLinear_t *SL, int i){
    int j, max = i;
    for(j = i+1; j < SL->n; j++){
               if(fabs(SL->A[j][i]) > fabs(SL->A[max][i]))
                    max = j;
    }
    return max;
}

//troca a linha de um sistema linear
void trocaLinha(SistLinear_t *SL, int i, int maior){
    int j = 0;
    real_t aux;

    //troca linhas da matriz
    for(j = 0; j < SL->n; j++){
        aux = SL->A[i][j];
        SL->A[i][j] = SL->A[maior][j];
        SL->A[maior][j] = aux;
    }

    //troca termos independentes
    aux = SL->b[i];
    SL->b[i] = SL->b[maior];
    SL->b[maior] = aux;
}

void retroS(SistLinear_t *SL, real_t *x){
    int i;
    for( i = SL->n - 1; i>= 0; i--){

        x[i] = SL->b[i];

        int j;
        for(j = i+1; j < SL->n; j++)
            x[i] -= SL->A[i][j] * x[j];
    
        x[i] /= SL->A[i][i];
    }   
}

void liberaVetor(real_t *x){
    //libera o vetor e nulifica os ponteiros
    free(x);
    x = NULL;
}

SistLinear_t *copiaMatriz(SistLinear_t *SL){
    SistLinear_t *SL2 = alocaSistLinear(SL->n);

    SL2->n = SL->n;
    SL2->erro = SL2->erro;

    int i, j;
    for(i = 0;i < SL->n;i++){
        for(j = 0;j < SL->n;j++)
            SL2->A[i][j] = SL->A[i][j];
    }
    
    for(i = 0;i < SL->n;i++)
        SL2->b[i] = SL->b[i];

    return SL2;
}

void pulaLinha(int n){
    int i;
    for(i = 0;i < n; i++)
        printf("\n");
}

real_t *residuo(SistLinear_t *SL, real_t *x){

    real_t *R = alocaVetor(SL->n);

    int i, j;
    for(i = 0; i < SL->n; i++){
        real_t AX = 0;
        for(j = 0; j < SL->n; j++){
            AX += SL->A[i][j] *x[j];
        }

        R[i] = SL->b[i] - AX;
    }

    return R;
}

void somaVetor(real_t *a, real_t *b, int n){
    int i;
    for( i = 0; i < n; i ++)
        a[i] += b[i];
}

void copiaVetor(real_t *a, real_t *b, int n){
    int i;
    for( i=0; i < n; i++)
        a[i] = b[i];
}

real_t maxDiff(real_t *a, real_t *b, int n){
    real_t max = b[0] - a[0];
    int i;
    for(i = 1; i < n; i++){
        if( max < (b[i] - a[i]) )
            max = b[i] - a[i];
    }

    return max;
}

void chamaRefinamento(SistLinear_t *SL, real_t *X, double *tTotal){
    SistLinear_t *SL3;
    real_t L2, *R;

    SL3 = copiaMatriz(SL);
    int k = refinamento(SL3, X, tTotal);

    printf("===> Refinamento: %f ms --> %d iteracoes\n--> X: ", *tTotal, k);
    prnVetor(X, SL->n);

    R = residuo(SL, X);
    L2 = normaL2Residuo(SL, X, R);
    printf("--> Norma L2 do residuo: %1.8e\n", L2);
    pulaLinha(1);

    liberaSistLinear(SL3);
}

int critDeConvergencia(SistLinear_t *SL){
    int i, j;
    real_t soma, max = 0.0;
    for(i = 0; i < SL->n; i++){

        soma = 0.0;
        //soma todos os elementos da linha
        for(j = 0; j < SL-> n; j++){
            if(i!=j)
                soma += fabs( SL->A[i][j] );
        }

        //divide pelo pivo
        soma /= fabs(SL->A[i][i]);
        if( soma > max)
            max = soma;
    }

    // alfa = max(alfai) < 1
    if( max > 1.0)
        return 0;
    return 1;
}
