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
*gauss(SYSTEM *A)//Функция, решающая СЛАУ A методом Гаусса.
{
    double temp;
    double **array = calloc(A->size, sizeof(*array));
    for (int i = 0; i < A->size; ++i) {
        array[i] = calloc(A->size, sizeof(**array));
    }
    double *vector = calloc(A->size, sizeof(*vector));
    for (int i = 0; i < A->size; ++i) {
        vector[i] = A->f[i];
        for (int j = 0; j < A->size; ++j) {
            array[i][j] = A->A[i][j];
        }
    }
    for (int j = 0; j < A->size; ++j) {//В данном цикле происходит приведение матрицы к верхней треугольной форме.
    	temp = array[j][j];
        vector[j] /= temp;
        for (int k = 0; k < A->size; ++k) {//Деление j-ой строки на j-ый диагональный элемент.
            array[j][k] /= temp;
        }
        for (int i = j + 1; i < A->size; ++i) {//Обнуление элементов j-ого столбца, стоящих ниже диагонали.
        	temp = array[i][j];
            vector[i] += -temp * vector[j];
            for (int k = 0; k < A->size; ++k) {
                array[i][k] += -temp * array[j][k];
            }
        }
    }
    double *result = calloc(A->size, sizeof(*result));
    for (int i = 0; i < A->size; ++i) {//В этом цикле происходит обратный ход (вычисление вектора решения).
        result[A->size - i - 1] = vector[A->size - i - 1];
        for (int j = A->size - i; j < A->size; ++j) {
            result[A->size - i - 1] -= array[A->size - i - 1][j];
        }
        for (int j = 0; j < A->size - i - 1; ++j) {
            array[j][A->size - i - 1] *= result[A->size - i - 1];
        }
    }
    for (int i = 0; i < A->size; ++i) {
        free(array[i]);
    }
    free(array);
    free(vector);
    return result;
}

double
*main_gauss(SYSTEM *A)//Функция, решающая СЛАУ A методом Гаусса с выбором главных компонент.
{
    double temp;
    double **array = calloc(A->size, sizeof(*array));
    for (int i = 0; i < A->size; ++i) {
        array[i] = calloc(A->size, sizeof(**array));
    }
    double *vector = calloc(A->size, sizeof(*vector));
    for (int i = 0; i < A->size; ++i) {
        vector[i] = A->f[i];
        for (int j = 0; j < A->size; ++j) {
            array[i][j] = A->A[i][j];
        }
    }
    int *columns = calloc(A->size, sizeof(*columns));
    for (int i = 0; i < A->size; ++i) {
        columns[i] = i;
    }
    for (int j = 0; j < A->size; ++j) {//В данном цикле происходит приведение матрицы к верхней треугольной форме.
        int i_max = j;
        for (int i = j; i < A->size; ++i) {//Поиск максимального элемента по j-ой строке верхней диагональной матрицы.
            if (fabs(array[j][i]) > fabs(array[j][i_max])) {
                i_max = i;
            }
        }
        if (i_max != j) {//Перестановка столбцов j и i_max.
            int tmp = columns[i_max];
            columns[i_max] = columns[j];
            columns[j] = tmp;
            for (int k = 0; k < A->size; ++k) {
                temp = array[k][j];
                array[k][j] = array[k][i_max];
                array[k][i_max] = temp;
            }
        }
        temp = array[j][j];
        vector[j] /= temp;
        for (int k = 0; k < A->size; ++k) {//Деление j-ой строки на её диагональный элемент.
            array[j][k] /= temp;
        }
        for (int i = j + 1; i < A->size; ++i) {//Обнуление элементов j-ого столбца, стоящих ниже диагонали.
        	temp = array[i][j];
            vector[i] += -temp * vector[j];
            for (int k = 0; k < A->size; ++k) {
                array[i][k] += -temp * array[j][k];
            }
        }
    }
    double *result = calloc(A->size, sizeof(*result));
    for (int i = 0; i < A->size; ++i) {//В этом цикле происходит обратный ход (вычисление вектора решения).
        result[A->size - i - 1] = vector[A->size - i - 1];
        for (int j = A->size - i; j < A->size; ++j) {
            result[A->size - i - 1] -= array[A->size - i - 1][j];
        }
        for (int j = 0; j < A->size - i - 1; ++j) {
            array[j][A->size - i - 1] *= result[A->size - i - 1];
        }
    }
    double *correct_result = calloc(A->size, sizeof(*correct_result));
    for (int i = 0; i < A->size; ++i) {
        correct_result[columns[i]] = result[i];
    }
    for (int i = 0; i < A->size; ++i) {
        free(array[i]);
    }
    free(columns);
    free(array);
    free(vector);
    free(result);
    return correct_result;
}

double
determinant(SYSTEM *A)//Функция, вычисляющая определитель матрицы A->A методом Гаусса с выбором главных компонент.
{
    double temp;
    double **array = calloc(A->size, sizeof(*array));
    for (int i = 0; i < A->size; ++i) {
        array[i] = calloc(A->size, sizeof(**array));
    }
    for (int i = 0; i < A->size; ++i) {
        for (int j = 0; j < A->size; ++j) {
            array[i][j] = A->A[i][j];
        }
    }
    double result = 1.0;//Определитель для матрицы array.
    for (int j = 0; j < A->size; ++j) {//В данном цикле происходит приведение матрицы к верхней треугольной форме.
        int i_max = j;
        for (int i = j; i < A->size; ++i) {//Поиск максимального элемента по j-ой строке верхней диагональной матрицы.
            if (fabs(array[j][i]) > fabs(array[j][i_max])) {
                i_max = i;
            }
        }
        if (i_max != j) {//Перестановка столбцов j и i_max.
            result *= -1;
            for (int k = 0; k < A->size; ++k) {
                temp = array[k][j];
                array[k][j] = array[k][i_max];
                array[k][i_max] = temp;
            }
        }
        temp = array[j][j];
        result *= temp;
        for (int k = 0; k < A->size; ++k) {//Деление j-ой строки на её диагональный элемент.
            array[j][k] /= temp;
        }
        for (int i = j + 1; i < A->size; ++i) {//Обнуление элементов j-ого столбца, стоящих ниже диагонали.
            temp = array[i][j];
            for (int k = 0; k < A->size; ++k) {
                array[i][k] += -temp * array[j][k];
            }
        }
    }
    for (int i = 0; i < A->size; ++i) {
        free(array[i]);
    }
    free(array);
    return result;
}

SYSTEM 
*reverse(SYSTEM *A)//Функция, находящая матрицу, обратную A->A, методом Гаусса-Жордана. 
{ 
    SYSTEM *result = calloc(1, sizeof(*result)); 
    result->size = A->size; 
    result->A = calloc(result->size, sizeof(*result->A));//Искомая матрица. 
    for (int i = 0; i < result->size; ++i) { 
        result->A[i] = calloc(result->size, sizeof(**result->A)); 
        for (int j = 0; j < result->size; ++j) { 
            result->A[i][j] = 0.0; 
        } 
        result->A[i][i] = 1.0; 
    } 
    double **array = calloc(A->size, sizeof(*array)); 
    for (int i = 0; i < A->size; ++i) { 
        array[i] = calloc(A->size, sizeof(**array)); 
        for (int j = 0; j < A->size; ++j) { 
            array[i][j] = A->A[i][j]; 
        } 
    } 
    for (int j = 0; j < A->size; ++j) {//В данном цикле происходит приведение матрицы array к единичной форме, рассматривая преобразования на матрице [A|I] => [I|revA]. 
        int i_max = j; 
        double temp; 
        for (int i = j; i < A->size; ++i) {//Поиск максимального элемента по j-му столбцу нижней диагональной матрицы, полученной из матрицы array. 
            if (fabs(array[i][j]) > fabs(array[i_max][j])) { 
                i_max = i; 
            } 
        } 
        if (i_max != j) {//Перестановка строк j и i_max. 
            for (int k = 0; k < A->size; ++k) { 
                temp = array[j][k]; 
                array[j][k] = array[i_max][k]; 
                array[i_max][k] = temp; 
                temp = result->A[j][k]; 
                result->A[j][k] = result->A[i_max][k]; 
                result->A[i_max][k] = temp; 
            } 
        } 
        temp = array[j][j]; 
        for (int k = 0; k < A->size; ++k) {//Деление новой j-ой строки на её диагональный элемент. 
            array[j][k] /= temp; 
            result->A[j][k] /= temp; 
        } 
        for (int i = 0; i < A->size; ++i) {//Обнуление элементов j-ого столбца, не являющихся диагональными. 
            if (i == j) { 
                continue; 
            } 
            temp = array[i][j]; 
            for (int k = 0; k < A->size; ++k) { 
                array[i][k] += -temp * array[j][k]; 
                result->A[i][k] += -temp * result->A[j][k]; 
            } 
        } 
    } 
    for (int i = 0; i < A->size; ++i) { 
        free(array[i]); 
    } 
    free(array); 
    return result; 
}

double
vect_norm(double *vector, int size)//Функция, вычисляющая евклидову норму вектора.
{
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        result += vector[i] * vector[i];
    }
    return sqrt(result);
}

double
frobenius_norm(double **matrix, int size)//Функция, вычисляющая фробениусовую норму матрицы.
{
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            result += matrix[i][j] * matrix[i][j];
        }
    }
    return sqrt(result);
}

double
cond_number(SYSTEM *A, SYSTEM *revA)//Функция, вычисляющая число обусловленности СЛАУ A.
{
    return frobenius_norm(A->A, A->size) * frobenius_norm(revA->A, revA->size);
}

void
analysis(double cond_A, double eps)//Оценка СЛАУ относительно числа обусловленности.
{
    if (cond_A > 100.0 || fabs(cond_A - 100.0) < eps) {
        puts("СЛАУ имеет плохую обусловленность");
    } else {
        puts("СЛАУ имеет хорошую обусловленность");
    }
}

double
rel_error(SYSTEM *A, double cond_A, double eps)//Верхняя оценка относительной погрешности решения СЛАУ A.
{
    double delta_vector = sqrt(A->size) * eps / vect_norm(A->f, A->size);
    double delta_matrix = A->size * eps / frobenius_norm(A->A, A->size);
    return (delta_vector + delta_matrix) * cond_A / (1 - cond_A * delta_matrix);
}

int
main(int argc, char **argv)
{
    double EPS;
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
        puts("Введите значение абсолютной погрешности вводимых данных");
        if (scanf("%lg", &EPS) <= 0) {
            perror("scanf()");
            return 1;
        }
        matrix->A = calloc(matrix->size, sizeof(*matrix->A));
        for (int i = 0; i < matrix->size; ++i) {
            matrix->A[i] = calloc(matrix->size, sizeof(**matrix->A));
        }
        matrix->f = calloc(matrix->size, sizeof(*matrix->f));
        puts("Введите (aij) и (bi)");
        for (int i = 0; i < matrix->size; ++i) {
            for (int j = 0; j < matrix->size; ++j) {
                if (scanf("%lg", &matrix->A[i][j]) <= 0) {
                    perror("scanf()");
                    return 1;
                }
            }
            if (scanf("%lg", &matrix->f[i]) <= 0) {
                perror("scanf()");
                return 1;
            }
        }
        matrix->x = gauss(matrix);
        puts("Решение СЛАУ методом Гаусса:");
        for (int i = 0; i < matrix->size; ++i) {
            printf("%.10g\n", matrix->x[i]);
        }
        free(matrix->x);
        matrix->x = main_gauss(matrix);
        puts("Решение СЛАУ методом Гаусса с выбором главных компонент:");
        for (int i = 0; i < matrix->size; ++i) {
            printf("%.10g\n", matrix->x[i]);
        }
        free(matrix->x);
        printf("Определитель: %.10g\n", determinant(matrix));
        SYSTEM *rev_matrix = reverse(matrix);
        puts("Обратная матрица:");
        for (int i = 0; i < rev_matrix->size; ++i) {
            for (int j = 0; j < rev_matrix->size; ++j) {
                printf("%.10g\t", rev_matrix->A[i][j]);
            }
            puts("");
        }
        double condition = cond_number(matrix, rev_matrix);
        printf("Число обусловленности СЛАУ: %.10g\n", condition);
        analysis(condition, EPS);
        printf("Относительная погрешность решения СЛАУ: ||dx||/||x|| <= %.10g\n", rel_error(matrix, condition, EPS));
        for (int i = 0; i < matrix->size; ++i) {
            free(matrix->A[i]);
            free(rev_matrix->A[i]);
        }
        free(matrix->A);
        free(rev_matrix->A);
        free(matrix->f);
        free(matrix);
        free(rev_matrix);
    }
    return 0;
}