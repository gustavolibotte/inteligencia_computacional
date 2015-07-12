#include <stdio.h>
#include <stdlib.h>

float* removeElementArray(float* oldArray, int arraySize, float removalElement)
{
    int i, j, ind;
    for(i = 0; i < arraySize; i++)
    {
        if(oldArray[i] == removalElement)
        {
            ind = i;
        }
    }
    float* newArray = realloc(oldArray, (arraySize - 1) * sizeof(float));
    j = 0;
    for(i = 0; i < arraySize; i++)
    {
        if(i != ind)
        {
            newArray[j] = oldArray[i];
            j++;
        }
    }
    return newArray;
}
