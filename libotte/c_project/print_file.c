#include <stdio.h>
#include <stdlib.h>
#include "print_file.h"

void printFile(float** finalPop, int bestInd, int dim, float fitness)
{
    FILE *f = fopen("C:\\Users\\Gustavo\\Documents\\Mestrado\\Inteligência Computacional\\Hammer\\Projeto\\Hammer\\resultados.txt", "a+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    int j;
    for(j = 0; j < dim; j++)
    {
        fprintf(f, "p[%d][%d] = %f\n", bestInd, j, finalPop[bestInd][j]);
    }
    fprintf(f, "\nfitness = %f\n", fitness);
    fprintf(f, "----------------------------------------------------------\n");

    fclose(f);
}
