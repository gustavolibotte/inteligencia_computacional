#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "normal.h"
#include "arrayManage.h"
#include "hammer.h"
#include "print_file.h"

//float f(float x[])
//{
//    return pow(pow(x[0], 2) + x[1] - 11, 2) + pow(x[0] + pow(x[1], 2) - 7, 2);
//}

struct Population
{
    float** weeds;
    float* fitness;
    int N_pop;
};

struct Population* createPopulation(float** weeds, int N_pop)
{
    struct Population* pop = malloc(sizeof(struct Population));
    pop->weeds = weeds;
    pop->N_pop = N_pop;
    return pop;
}

struct Population* initPop(int N_0, int dim, float Lb[][dim], float Ub[][dim])
{
    int i, j;
    float** pop = malloc(N_0 * sizeof(int*));
    for(i = 0; i < N_0; i++)
    {
        pop[i] = malloc(dim * sizeof(int));
    }
    srand(time(NULL));
    for(i = 0; i < N_0; i++)
    {
        for(j = 0; j < dim; j++)
        {
            if(j == 6)
            {
                float tempFuelMat = Lb[0][j] + rand() / (RAND_MAX / (Ub[0][j] - Lb[0][j]) + 1);
                if(tempFuelMat < 100.5)
                {
                    pop[i][j] = 100;
                }
                else
                {
                    pop[i][j] = 101;
                }
            }
            else if(j == 7)
            {
                float tempCladdingMat = Lb[0][j] + rand() / (RAND_MAX / (Ub[0][j] - Lb[0][j]) + 1);
                if(tempCladdingMat < 200.5)
                {
                    pop[i][j] = 200;
                }
                else if(tempCladdingMat >= 200.5 && tempCladdingMat < 201.5)
                {
                    pop[i][j] = 201;
                }
                else
                {
                    pop[i][j] = 202;
                }
            }
            else
            {
                float tempParam = Lb[0][j] + rand() / (RAND_MAX / (Ub[0][j] - Lb[0][j]) + 1);
                while(tempParam < Lb[0][j] | tempParam > Ub[0][j])
                {
                    tempParam = Lb[0][j] + rand() / (RAND_MAX / (Ub[0][j] - Lb[0][j]) + 1);
                }
                pop[i][j] = tempParam;
            }
        }
    }
    struct Population *newPop = createPopulation(pop, N_0);
    return newPop;
}

float getMinMaxFitness(struct Population* pop, int optFlag)
{
    int i;
    int N_pop = pop->N_pop;
    float* fitness = pop->fitness;
    float ansFitness = fitness[0];
    for(i = 0; i < N_pop; i++)
    {
        if(optFlag == 0)
        {
            if(fitness[i] < ansFitness)
            {
                ansFitness = fitness[i];
            }
        }
        else if(optFlag == 1)
        {
            if(fitness[i] > ansFitness)
            {
                ansFitness = fitness[i];
            }
        }
    }
    return ansFitness;
}

int* reproduction(struct Population* pop, int s_min, int s_max, int dim)
{
    int i;
    int numWeeds = pop->N_pop;
    float** currentWeeds = pop->weeds;
    float* currentFitness = malloc(numWeeds * sizeof(int));
    for(i = 0; i < numWeeds; i++)
    {
        currentFitness[i] = eval(currentWeeds[i]);
    }
    pop->fitness = currentFitness;
    float oppositeKathete = s_max - s_min;
    float adjacentKathete = getMinMaxFitness(pop, 1) - getMinMaxFitness(pop, 0);
    float angularCoefficient = atan(oppositeKathete / adjacentKathete);
    float linearCoefficient = s_min - (angularCoefficient * getMinMaxFitness(pop, 0));
    int* seeds = malloc(numWeeds * sizeof(int));
    for(i = 0; i < numWeeds; i++)
    {
        seeds[i] = (int)(s_max - ((angularCoefficient * eval(currentWeeds[i])) + linearCoefficient));
    }
    return seeds;
}

struct Population* spatialDispersal(struct Population* pop, int dim, float sigma_initial, float sigma_final, float n, int it_max, int it, int* seeds, float Lb[][dim], float Ub[][dim])
{
    int i, j, k, iterator;
    int randomSeedGenerator = 123456789;
    int numWeeds = pop->N_pop;
    float** currentWeeds = pop->weeds;
    float mu = 0;
    float sigma_it = ((pow(it_max - it, n) / pow(it_max, n)) * (sigma_initial - sigma_final)) + sigma_final;
    int totalSeeds = 0;
    for(i = 0; i < numWeeds; i++)
    {
        totalSeeds += seeds[i];
    }
    float** newPlants = malloc(totalSeeds * sizeof(int*));
    for(i = 0; i < totalSeeds; i++)
    {
        newPlants[i] = malloc(dim * sizeof(int));
    }
    iterator = 0;
    float tempNewWeed = 0;
    for(i = 0; i < numWeeds; i++)
    {
        if(seeds[i] > 0)
        {
            for(j = 0; j < seeds[i]; j++)
            {
                for(k = 0; k < dim; k++)
                {
                    float r = r8_normal_ab(mu, sigma_it, &randomSeedGenerator);
                    if(((float)rand() / (float)RAND_MAX) > 0.5)
                    {
                        tempNewWeed = currentWeeds[i][k] + r;
                    }
                    else
                    {
                        tempNewWeed = currentWeeds[i][k] - r;
                    }

                    if(k >= 0 && k <= 5)
                    {
                        if(tempNewWeed >= Lb[0][k] && tempNewWeed <= Ub[0][k])
                        {
                            newPlants[iterator][k] = tempNewWeed;
                        }
                        else
                        {
                            newPlants[iterator][k] = currentWeeds[i][k];
                        }
                    }
                    else
                    {
                        if(k == 6)
                        {
                            if(tempNewWeed < 100.5)
                            {
                                newPlants[iterator][k] = 100;
                            }
                            else
                            {
                                newPlants[iterator][k] = 101;
                            }
                        }
                        else
                        {
                            if(tempNewWeed < 200.5)
                            {
                                newPlants[iterator][k] = 200;
                            }
                            else if(tempNewWeed >= 200.5 && tempNewWeed < 201.5)
                            {
                                newPlants[iterator][k] = 201;
                            }
                            else
                            {
                                newPlants[iterator][k] = 202;
                            }
                        }
                    }
                }
                iterator++;
            }
        }
    }
    struct Population *newPop = createPopulation(newPlants, totalSeeds);
    return newPop;
}

struct Population* competitiveExclusion(struct Population* pop, struct Population* newPlants, int dim, int p_max)
{
    int i, j, k, w, ind;
    int numPop = pop->N_pop + newPlants->N_pop;
    float minTempFitness;
    int* removalIndex;
    float** oldPopulation = pop->weeds;
    float** newPopulation = newPlants->weeds;
    float** tempPopulation = malloc(numPop * sizeof(int*));
    for(i = 0; i < numPop; i++)
    {
        tempPopulation[i] = malloc(dim * sizeof(int));
    }
    j = 0;
    for(i = 0; i < pop->N_pop; i++)
    {
        for(k = 0; k < dim; k++)
        {
            tempPopulation[j][k] = oldPopulation[i][k];
        }
        j++;
    }
    for(i = 0; i < newPlants->N_pop; i++)
    {
        for(k = 0; k < dim; k++)
        {
            tempPopulation[j][k] = newPopulation[i][k];
        }
        j++;
    }
    float* fitnessTempPopulation = malloc(numPop * sizeof(int));
    for(i = 0; i < numPop; i++)
    {
        fitnessTempPopulation[i] = eval(tempPopulation[i]);
    }
    float** finalPopulation = malloc(p_max * sizeof(int*));
    for(i = 0; i < p_max; i++)
    {
        finalPopulation[i] = malloc(dim * sizeof(int));
    }
    for(i = 0; i < p_max; i++)
    {
        minTempFitness = fitnessTempPopulation[0];
        for(j = 0; j < (numPop - i); j++)
        {
            if(fitnessTempPopulation[j] <= minTempFitness)
            {
                minTempFitness = fitnessTempPopulation[j];
                ind = j;
            }
        }
        if(i == 0)
        {
            removalIndex = malloc(1 * sizeof(int));
            removalIndex[0] = ind;
        }
        else
        {
            removalIndex = realloc(removalIndex, (i + 1) * sizeof(int));
            removalIndex[i] = ind;
        }
        fitnessTempPopulation = removeElementArray(fitnessTempPopulation, (numPop - i), fitnessTempPopulation[ind]);
    }
    for(i = 0; i < p_max; i++)
    {
        finalPopulation[i] = tempPopulation[removalIndex[i] + i];
    }
    return createPopulation(finalPopulation, p_max);
}

void iwo()
{
    int it, i, j;
    int dim = 8;
    int N_0 = 10;
    int p_max = 15;
    int s_min = 0;
    int s_max = 5;
    int it_max = 20;
    float sigma_initial = 1;
    float sigma_final = 0.01;
    float n = 3;
    float Lb[][8] = {0.508, 0.025, 0.025, 2, 2, 2, 100, 200};
    float Ub[][8] = {1.27, 0.254, 0.762, 5, 5, 5, 101, 202};

    int numRuns = 100;
    int iRuns;
    for(iRuns = 0; iRuns < numRuns; iRuns++)
    {
        struct Population* pop = initPop(N_0, dim, Lb, Ub);
        for(it = 0; it < it_max; it++)
        {
            int* seeds = reproduction(pop, s_min, s_max, dim);
            struct Population* newPlants = spatialDispersal(pop, dim, sigma_initial, sigma_final, n, it_max, it, seeds, Lb, Ub);
            pop = competitiveExclusion(pop, newPlants, dim, p_max);
        }
        float** finalPopulation = pop->weeds;

        int bestIndex = 0;
        float bestFitness = eval(finalPopulation[0]);
        for(i = 0; i < pop->N_pop; i++)
        {
            if(eval(finalPopulation[i]) <= bestFitness)
            {
                bestFitness = eval(finalPopulation[i]);
                bestIndex = i;
            }

        }
        printf("Iteracao = %d", iRuns);
        printFile(finalPopulation, bestIndex, dim, bestFitness);
    }
}

int main()
{
    iwo();
    return 0;
}
