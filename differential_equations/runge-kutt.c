#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct net
{
    int ynum;
    int count;
    double h;
    double **net;//net[count][ynum + 1]; значение: 0 - x, i - yi
} NET;

typedef struct system
{
    int n;
    double (**f)(NET *, double *, int);//f[n], double * - альтернативные аргументы, int - текущая итерация
} SYSTEM;

double
f1(NET *x, double *vect, int i)//SYSTEM 2-21 (x, y1, y2) = (0, 1, 0.25)
{
    double result = 0;
    if (x != NULL) {//2.4 * v - u
        result = 2.4 * x->net[i][2] - x->net[i][1];
    } else {
        result = 2.4 * vect[2] - vect[1];
    }
    return result;
}

double
f2(NET *x, double *vect, int i)//SYSTEM 2-21
{
    double result = 0;
    if (x != NULL) {//e^(-u) - x + 2.2 * v
        result = exp(-x->net[i][1]) - x->net[i][0] + 2.2 * x->net[i][2];
    } else {
        result = exp(-vect[1]) - vect[0] + 2.2 * vect[2];
    }
    return result;
}

void
runge_kutt_2(NET *args, SYSTEM *funcs, double a)
{
    double *tmp = calloc(args->ynum + 1, sizeof(*tmp));
    for (int i = 1; i < args->count; ++i) {
        tmp[0] = args->net[i - 1][0] + args->h/(2 * a);
        for (int j = 1; j < args->ynum + 1; ++j) {
            tmp[j] = args->net[i - 1][j] + (args->h/(2 * a)) * funcs->f[j - 1](args, NULL, i - 1);
        }
        for (int j = 1; j < args->ynum + 1; ++j) {
            args->net[i][j] = args->net[i - 1][j] + args->h * ((1 - a) * funcs->f[j - 1](args, NULL, i - 1) + a * funcs->f[j - 1](NULL, tmp, i - 1));
        }
    }
    free(tmp);
    return;
}

void
runge_kutt_4(NET *args, SYSTEM *funcs)
{
    double **tmp = calloc(3, sizeof(*tmp));//tmp[3][ynum + 1] 0: аргументы для p2, 1: аргументы для p3, 2: аргументы для p4
    double **p = calloc(3, sizeof(*p));
    for (int i = 0; i < 3; ++i) {
        p[i] = calloc(args->ynum, sizeof(**p));
    }
    for (int i = 0; i < 3; ++i) {
        tmp[i] = calloc(args->ynum + 1, sizeof(**tmp));
    }
    for (int i = 1; i < args->count; ++i) {
        tmp[0][0] = args->net[i - 1][0] + args->h/2;
        for (int j = 1; j < args->ynum + 1; ++j) {
            tmp[0][j] = args->net[i - 1][j] + (args->h/2) * funcs->f[j - 1](args, NULL, i - 1);
        }//Аргументы для p2 вычислены
        for (int j = 0; j < args->ynum; ++j) {
            p[0][j] = funcs->f[j](NULL, tmp[0], i - 1);
        }//p2 вычислены
        tmp[1][0] = tmp[0][0];
        for (int j = 1; j < args->ynum + 1; ++j) {
            tmp[1][j] = args->net[i - 1][j] + (args->h/2) * p[0][j - 1];
        }//Аргументы для p3 вычислены
        for (int j = 0; j < args->ynum; ++j) {
            p[1][j] = funcs->f[j](NULL, tmp[1], i - 1);
        }//p3 вычислены
        tmp[2][0] = args->net[i - 1][0] + args->h;
        for (int j = 1; j < args->ynum + 1; ++j) {
            tmp[2][j] = args->net[i - 1][j] + args->h * p[1][j - 1];
        }//Аргументы для p4 вычислены
        for (int j = 0; j < args->ynum; ++j) {
            p[2][j] = funcs->f[j](NULL, tmp[2], i - 1);
        }//p4 вычислены
        for (int j = 1; j < args->ynum + 1; ++j) {
            args->net[i][j] = args->net[i - 1][j] + (args->h/6) * (funcs->f[j - 1](args, NULL, i - 1) + 2 * p[0][j - 1] + 2 * p[1][j - 1] + p[2][j - 1]);
        }//y i-ой итерации вычислены
    }
    for (int i = 0; i < 3; ++i) {
        free(p[i]);
    }
    free(p);
    for (int i = 0; i < 3; ++i) {
        free(tmp[i]);
    }
    free(tmp);
    return;
}

int
main(int argc, char **argv)
{
    int n;
    puts("Введите число уравнений в системе:");
    if (scanf("%d", &n) < 1) {
        perror("scanf()");
        return 1;
    }
    if (n > 2) {
        perror("Превышено число уравнений в системе");
        return 2;
    }
    double h;
    puts("Введите шаг итерации:");
    if (scanf("%lf", &h) < 1) {
        perror("scanf()");
        return 1;
    }
    double x0, b;
    puts("Введите начальные условия (x0, b, y01, ..., y0n):");
    if (scanf("%lf%lf", &x0, &b) < 2) {
        perror("scanf()");
        return 1;
    }
    NET *solution = calloc(1, sizeof(*solution));
    solution->ynum = n;
    int count = (b - x0)/h + 1;
    solution->count = count;
    solution->h = h;
    solution->net = calloc(count, sizeof(*solution->net));
    for (int i = 0; i < count; ++i) {
        solution->net[i] = calloc(n + 1, sizeof(**solution->net));
    }
    solution->net[0][0] = x0;
    for (int i = 1; i < count - 1; ++i) {
        solution->net[i][0] = solution->net[i - 1][0] + h;
    }
    solution->net[count - 1][0] = b;
    for (int i = 0; i < n; ++i) {
        if (scanf("%lf", &solution->net[0][i + 1]) < 1) {
            perror("scanf()");
            return 1;
        }
    }
    SYSTEM *equations = calloc(1, sizeof(*equations));
    equations->n = n;
    equations->f = calloc(2, sizeof(*equations->f));
    equations->f[0] = f1;
    equations->f[1] = f2;
    runge_kutt_2(solution, equations, 0.5);
    puts("Решение методом Рунге-Кутта 2 порядка точности с параметром 0.5");
    for (int j = 1; j < n + 1; ++j) {
        printf("ListPlot[{");
        for (int i = 0; i < count - 1; ++i) {
            printf("{%.10g, %.10g}, ", solution->net[i][0], solution->net[i][j]);
        }
        printf("{%.10g, %.10g}}]\n", solution->net[count - 1][0], solution->net[count - 1][j]);
    }
    runge_kutt_2(solution, equations, 1);
    puts("Решение методом Рунге-Кутта 2 порядка точности с параметром 1");
    for (int j = 1; j < n + 1; ++j) {
        printf("ListPlot[{");
        for (int i = 0; i < count - 1; ++i) {
            printf("{%.10g, %.10g}, ", solution->net[i][0], solution->net[i][j]);
        }
        printf("{%.10g, %.10g}}]\n", solution->net[count - 1][0], solution->net[count - 1][j]);
    }
    runge_kutt_4(solution, equations);
    puts("Решение методом Рунге-Кутта 4 порядка точности");
    for (int j = 1; j < n + 1; ++j) {
        printf("ListPlot[{");
        for (int i = 0; i < count - 1; ++i) {
            printf("{%.10g, %.10g}, ", solution->net[i][0], solution->net[i][j]);
        }
        printf("{%.10g, %.10g}}]\n", solution->net[count - 1][0], solution->net[count - 1][j]);
    }
    free(equations->f);
    free(equations);
    for (int i = 0; i < count; ++i) {
        free(solution->net[i]);
    }
    free(solution);
    return 0;
}