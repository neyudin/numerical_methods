#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct system//Структура, хранящая СЛАУ.
{
    int size;
    double **A;
    double *x;
    double *f;
} SYSTEM;

double
vect_norm(double *vect, int size)//Функция, вычисляющая евклидову норму вектора.
{
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        result += vect[i] * vect[i];
    }
    return sqrt(result);
}

int
exact_converge(double *prev, SYSTEM *A, double *eps)//Функция, вычисляющая сходимость решения по невязке левой и правой частей СЛАУ.
{
    double *dx = calloc(A->size, sizeof(*dx));
    for (int i = 0; i < A->size; ++i) {
        dx[i] = prev[i] - A->x[i];
    }
    double result = vect_norm(dx, A->size);
    free(dx);
    if (result < *eps) {
        return 1;
    }
    return 0;
}

unsigned long long
seidel(SYSTEM *A, double *eps)//Функция, решающая СЛАУ итерационным методом Зейделя.
{
    double *prev = calloc(A->size, sizeof(*prev));
    unsigned long long counter = 0;
    do {
        for (int i = 0; i < A->size; ++i) {
            prev[i] = A->x[i];
        }
        for (int i = 0; i < A->size; ++i) {
            double temp = 0.0;
            for (int j = 0; j < i; ++j) {
                temp += A->A[i][j] * A->x[j];
            }
            for (int j = i + 1; j < A->size; ++j) {
                temp += A->A[i][j] * prev[j];
            }
            A->x[i] = (A->f[i] - temp) / A->A[i][i];
        }
        ++counter;
    } while (!exact_converge(prev, A, eps));
    free(prev);
    return counter;
}

unsigned long long
upper_relaxation(SYSTEM *A, double *omega, double *eps)//Функция, решающая СЛАУ итерационным методом верхней релаксации.
{
    double *prev = calloc(A->size, sizeof(*prev));
    unsigned long long counter = 0;
    do {
        for (int i = 0; i < A->size; ++i) {
            prev[i] = A->x[i];
        }
        for (int i = 0; i < A->size; ++i) {
            double temp = 0.0;
            for (int j = 0; j < i; ++j) {
                temp += A->A[i][j] * A->x[j];
            }
            for (int j = i; j < A->size; ++j) {
                temp += A->A[i][j] * prev[j];
            }
            A->x[i] = prev[i] + *omega * (A->f[i] - temp) / A->A[i][i];
        }
        ++counter;
    } while (!exact_converge(prev, A, eps));
    free(prev);
    return counter;
}

int
main(int argc, char **argv)
{
    double EPS, OMEGA;
    int count;
    puts("Введите количество СЛАУ, которых необходимо решить");
    if (scanf("%d", &count) <= 0) {
        perror("scanf()");
        return 1;
    }
    for (int k = 0; k < count; ++k) {
        printf("\n==========\t\t\tСЛАУ номер %d\t\t\t==========\n\n", k + 1);
        SYSTEM *matrix = calloc(1, sizeof(*matrix));
        puts("Введите размер СЛАУ (n)");
        if (scanf("%d", &matrix->size) <= 0) {
            perror("scanf()");
            return 1;
        }
        puts("Введите значение итерационного параметра для метода верхней релаксации");
        if (scanf("%lg", &OMEGA) <= 0) {
            perror("scanf()");
            return 1;
        }
        puts("Введите значение абсолютной погрешности вводимых данных");
        if (scanf("%lg", &EPS) <= 0) {
            perror("scanf()");
            return 1;
        }
        matrix->A = calloc(matrix->size, sizeof(*matrix->A));
        for (int i = 0; i < matrix->size; ++i) {
            matrix->A[i] = calloc(matrix->size, sizeof(**matrix->A));
        }
        matrix->x = calloc(matrix->size, sizeof(*matrix->x));
        matrix->f = calloc(matrix->size, sizeof(*matrix->f));
        puts("Введите (aij) и (bi)");
        for (int i = 0; i < matrix->size; ++i) {
            for (int j = 0; j < matrix->size; ++j) {
                if (scanf("%lg", &matrix->A[i][j]) <= 0) {
                    perror("scanf()");
                    return 1;
                }
            }
            matrix->x[i] = 0.0;
            if (scanf("%lg", &matrix->f[i]) <= 0) {
                perror("scanf()");
                return 1;
            }
        }
        unsigned long long count_1 = seidel(matrix, &EPS);
        puts("Решение СЛАУ методом Зейделя:");
        for (int i = 0; i < matrix->size; ++i) {
            printf("%.10g\n", matrix->x[i]);
        }
        printf("Было выполнено %llu итераций\n", count_1);
        for (int i = 0; i < matrix->size; ++i) {
            matrix->x[i] = 0.0;
        }
        unsigned long long count_2 = upper_relaxation(matrix, &OMEGA, &EPS);
        puts("Решение СЛАУ методом верхней релаксации:");
        for (int i = 0; i < matrix->size; ++i) {
            printf("%.10g\n", matrix->x[i]);
        }
        printf("Было выполнено %llu итераций\n", count_2);
        for (int i = 0; i < matrix->size; ++i) {
            free(matrix->A[i]);
        }
        free(matrix->A);
        free(matrix->x);
        free(matrix->f);
        free(matrix);
    }
    return 0;
}
