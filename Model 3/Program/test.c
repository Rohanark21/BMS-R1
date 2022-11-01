#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

int main(void)
{
    int matdata=fopen("PANdata_P25.mat","r");
    if (matdata == NULL)
    {
        fprintf(stderr, "Unable to open file %s\n", "PANdata_P25.mat");
        return 1; /* Exit to operating system */
    }
    return 0;
}


