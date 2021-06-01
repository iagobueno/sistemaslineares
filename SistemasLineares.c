#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "functions.h"

/*!
\brief Essa função calcula a norma L2 do resíduo de um sistema linear 

\param SL Ponteiro para o sistema linear
\param x Solução do sistema linear
\param res Valor do resíduo

\return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
}

/*!
\brief Método da Eliminação de Gauss

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução
\param tTotal tempo gasto pelo método

\return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal)
{
    //percorre as linhas da matriz
    int i;
    for(i = 0; i < (SL->n); i++){

        //pivoteamento parcial
        int maior = encontraMaior(SL, i);
        if(i != maior)
            trocaLinha(SL, i, maior);

        int k;
        for(k = i+1; k < SL->n; k++){
            
            //calcula M e zera o elemento
            double m = SL->A[k][i] / SL->A[i][i];
            SL->A[k][i] = 0;

            //percorre as colunas da matriz
            int j;
            for(j = i+1; j < SL->n; j++){
                SL->A[k][j] -= SL->A[i][j] * m;
            }

            //termo independente
            SL->b[k] -= SL->b[i] * m;
        }
    }
    retroS(SL, x);
}

/*!
\brief Método de Jacobi

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int gaussJacobi(SistLinear_t *SL, real_t *x, double *tTotal)
{
}

/*!
\brief Método de Gauss-Seidel

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int gaussSeidel(SistLinear_t *SL, real_t *x, double *tTotal)
{
}

/*!
\brief Método de Refinamento

\param SL Ponteiro para o sistema linear
\param x ponteiro para o vetor solução. Ao iniciar função contém
    valor inicial para início do refinamento
\param tTotal tempo gasto pelo método

\return código de erro. Um nr positivo indica sucesso e o nr
de iterações realizadas. Um nr. negativo indica um erro:
-1 (não converge) -2 (sem solução)
*/
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal)
{
}

/*!
\brief Alocaçao de memória 

\param n tamanho do SL

\return ponteiro para SL. NULL se houve erro de alocação
*/
SistLinear_t *alocaSistLinear(unsigned int n)
{
    SistLinear_t *SL = malloc ( sizeof(SistLinear_t));
    testaMalloc(SL);

    SL->n = n;

    SL->A = alocaCoef(n);
    SL->b = alocaVetor(n);

    return SL;
}

/*!
\brief Liberaçao de memória 

\param sistema linear SL
*/
void liberaSistLinear(SistLinear_t *SL)
{
    free(SL->b);

    free(SL->A[0]);
    free(SL->A);
    free(SL);
}

/*!
\brief Leitura de SL a partir de Entrada padrão (stdin).

\return sistema linear SL. NULL se houve erro (leitura ou alocação)
*/
SistLinear_t *lerSistLinear()
{
    int n, fim;
    real_t erro;

    fim = scanf("%d%f", &n, &erro);
    if(fim == EOF)
        exit(0);

    SistLinear_t *SL = alocaSistLinear(n);

    SL->erro = erro;

    int i, j;
    for(i = 0;i < n;i++){
        for(j = 0;j < n;j++)
            scanf("%f", SL->A[i]+j);
    }
    
    for(i = 0;i < n;i++)
        scanf("%f", SL->b+i);

    return SL;
}

// Exibe SL na saída padrão
void prnSistLinear(SistLinear_t *SL)
{
    int i, j;
    for(i = 0; i < SL->n; i++){
        for(j = 0; j < SL->n; j++)
            printf("%f ", SL->A[i][j]); 
        printf("\n");   
    }
}

// Exibe um vetor na saída padrão
void prnVetor(real_t *v, unsigned int n)
{
    int i;
    for(i = 0; i < n; i++){
        printf("%f ", v[i]);
    }
    printf("\n");
}
