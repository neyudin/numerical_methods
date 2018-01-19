#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void
create_test(FILE *f, int N, int x, int y, int z)//Функция, создающая СЛАУ для итерационных методов.
{
    double **array = calloc(N, sizeof(*array));
    double *vector = calloc(N, sizeof(*vector));
    for (int i = 0; i < N; ++i) {
        array[i] = calloc(N, sizeof(**array));
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            array[i][j] = rand() / (RAND_MAX + 1.0) * x;
            array[j][i] = array[i][j];
        }
        vector[i] = rand() / (RAND_MAX + 1.0) * y - y / 2.0;
    }
    for (int i = 0; i < N; ++i) {
        double temp = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                continue;
            }
            temp += array[i][j];
        }
        array[i][i] = temp + rand() / (RAND_MAX + 1.0) * z + 1.0;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, "%.10g ", array[i][j]);
        }
        fprintf(f, "%.10g\n", vector[i]);
    }
    for (int i = 0; i < N; ++i) {
        free(array[i]);
    }
    free(array);
    free(vector);
}

int
main(int argc, char **argv)//Данный генератор тестов создаёт тест из одой СЛАУ для итерационных методов, записанной в отдельном файле.
{
	/*Аргументы командной строки:
	argv[1] - Размер СЛАУ.
	argv[2] - Итерационный параметр.
	argv[3] - Абсолютная погрешность вводимых данных.
	argv[4] - "Семя" для генератора псевдослучайных чисел.
    argv[5] - Параметр x для функции, создающей СЛАУ.
    argv[6] - Параметр y для функции, создающей СЛАУ.
    argv[7] - Параметр z для функции, создающей СЛАУ.
	argv[8] - Имя файла, в который будет записана СЛАУ.
	*/
    if (argc < 9) {
        perror("Недостачно аргументов командной строки");
        return 1;
    }
    FILE *file = fopen(argv[8], "a");
    if (!file) {
        perror("fopen()");
        return 2;
    }
    int size = strtol(argv[1], NULL, 10);
    fprintf(file, "%d\n%d\n%.10g\n%.10g\n", 1, size, strtod(argv[2], NULL), strtod(argv[3], NULL));
    srand(strtol(argv[4], NULL, 10));
    create_test(file, size, strtol(argv[5], NULL, 10), strtol(argv[6], NULL, 10), strtol(argv[7], NULL, 10));
    if (fclose(file)) {
        perror("fclose()");
        return 2;
    }
    return 0;
}
