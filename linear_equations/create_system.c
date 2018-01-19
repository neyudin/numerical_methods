#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum
{
    M = 6
};

typedef struct system//Структура, хранящая СЛАУ.
{
    int size;
    double **A;
    double *f;
} SYSTEM;

SYSTEM
*create_system(int size, double x)//Данная функция создаёт СЛАУ.
{
    SYSTEM *result = calloc(1, sizeof(*result));
    result->size = size;
    result->A = calloc(result->size, sizeof(*result->A));
    result->f = calloc(result->size, sizeof(*result->f));
    for (int i = 0; i < result->size; ++i) {
        result->A[i] = calloc(result->size, sizeof(**result->A));
    }
    double qm = 1.001 - 2 * M * 0.001;
    for (int i = 1; i <= result->size; ++i) {
        result->f[i - 1] = x * exp(x/i) * cos(x/i);
        for (int j = 1; j <= result->size; ++j) {
            if (i == j) {
                result->A[i - 1][j - 1] = pow(qm - 1.0, (double) (i + j));
            } else {
                result->A[i - 1][j - 1] = pow(qm, (double) (i + j)) + 0.1 * (j - i);
            }
        }
    }
    return result;
}

int
main(int argc, char **argv)//Данный генератор тестов создаёт тест из одой СЛАУ, записанной в отдельном файле.
{
    /*Аргументы командной строки:
    argv[1] - Размер СЛАУ.
    argv[2] - Абсолютная пргрешность вводимых данных.
    argv[3] - Параметр x для функции, создающей СЛАУ.
    argv[4] - Имя файла, в который будет записана СЛАУ.
    */
    if (argc < 5) {
        perror("Недостаточно аргументов командной строки");
        return 1;
    }
    FILE *f = fopen(argv[4], "a");
    if (!f) {
        perror("fopen()");
        return 2;
    }
    SYSTEM *A = create_system(strtol(argv[1], NULL, 10), strtod(argv[3], NULL));
    fprintf(f, "%d\n%d\n%.10g\n", 1, A->size, strtod(argv[2], NULL));
    for (int i = 0; i < A->size; ++i) {
        for (int j = 0; j < A->size; ++j) {
            fprintf(f, "%.10g ", A->A[i][j]);
        }
        fprintf(f, "%.10g\n", A->f[i]);
    }
    if (fclose(f)) {
        perror("fclose()");
        return 2;
    }
    for (int i = 0; i < A->size; ++i) {
        free(A->A[i]);
    }
    free(A->A);
    free(A->f);
    free(A);
    return 0;
}