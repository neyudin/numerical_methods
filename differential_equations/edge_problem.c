#include <stdio.h>
#include <stdlib.h>
#include <math.h>

volatile double EPS = 0.00000000001;

typedef struct system
{
    double (*p)(double);
    double (*q)(double);
    double (*f)(double);
} SYSTEM;

typedef struct condition
{
    double bounds[2];
    double sigma[2];
    double gamma[2];
    double delta[2];
} COND;

typedef struct net
{
    double h;
    int count;
    double **net;
} NET;

void
edge_problem(SYSTEM *equation, COND *values, NET *args)
{
    double **coefs = calloc(args->count - 1, sizeof(*coefs));
    for (int i = 0; i < args->count - 1; ++i) {
        coefs[i] = calloc(2, sizeof(**coefs));
    }
    coefs[0][0] = values->gamma[0]/(values->gamma[0] - args->h * values->sigma[0]);
    coefs[0][1] = args->h * values->delta[0]/(args->h * values->sigma[0] - values->gamma[0]);
    for (int i = 1; i < args->count - 1; ++i) {
        coefs[i][0] = -(2 + equation->p(args->net[i - 1][0]) * args->h)/((2 - equation->p(args->net[i - 1][0]) * args->h) * coefs[i - 1][0] - 4 + 2 * equation->q(args->net[i - 1][0]) * args->h * args->h);
        coefs[i][1] = -(2 * equation->f(args->net[i - 1][0]) * args->h * args->h + (2 - equation->p(args->net[i - 1][0]) * args->h) * coefs[i - 1][1])/((2 - equation->p(args->net[i - 1][0]) * args->h) * coefs[i - 1][0] - 4 + 2 * equation->q(args->net[i - 1][0]) * args->h * args->h);
    }
    if (fabs(values->gamma[1]) < EPS) {
        args->net[args->count - 1][1] = values->delta[1]/values->sigma[1];
    } else {
        if (fabs(values->sigma[1]) < EPS) {
            args->net[args->count - 1][1] = (coefs[args->count - 2][1] + values->delta[1] * args->h/values->gamma[1])/(1 - coefs[args->count - 2][0]);
        } else {
            args->net[args->count - 1][1] = (coefs[args->count - 2][1] + values->delta[1] * args->h/values->gamma[1])/(values->sigma[1] * args->h/values->gamma[1] + 1- coefs[args->count - 2][0]);
        }
    }
    for (int i = 1; i < args->count; ++i) {
        args->net[args->count - i - 1][1] = coefs[args->count - i - 1][0] * args->net[args->count - i][1] + coefs[args->count - i - 1][1];
    }
    for (int i = 0; i < args->count - 1; ++i){
        free(coefs[i]);
    }
    free(coefs);
    return;
}

double
p(double x)
{
    double result = 2.0;//4
    //double result = -1.0;//6
    //double result = -0.5 * x;//12
    return result;
}

double
q(double x)
{
    double result = -1.0/x;//4
    //double result = 2.0/x;//6
    //double result = 1.0;//12
    return result;
}

double
f(double x)
{
    double result = -3.0;//4
    //double result = -x - 0.4;//6
    //double result = -2.0;//12
    return result;
}

int
main(int argc, char **argv)
{
    SYSTEM equation;
    equation.p = p;
    equation.q = q;
    equation.f = f;
    COND values;
    puts("Введите условия краевой задачи для ОДУ 2 порядка: y\" + p(x) * y' + q(x) * y = -f(x)");
    puts("s1 * y(a) + g1 * y'(a) = d1");
    puts("s2 * y(b) + g2 * y'(b) = d2");
    puts("В формате:a,b,s1,g1,d1,s2,g2,d2");
    if (scanf("%lf%lf%lf%lf%lf%lf%lf%lf", &values.bounds[0], &values.bounds[1], &values.sigma[0], &values.gamma[0], &values.delta[0], &values.sigma[1], &values.gamma[1], &values.delta[1]) < 8) {
        perror("scanf()");
        return 1;
    }
    NET *args = calloc(1, sizeof(*args));
    puts("Введите шаг итерации(>10^(-11))");
    if (scanf("%lf", &args->h) < 1) {
        perror("scanf()");
        return 1;
    }
    args->count = (values.bounds[1] - values.bounds[0])/args->h + 1;
    args->net = calloc(args->count, sizeof(*args->net));
    for (int i = 0; i < args->count; ++i) {
        args->net[i] = calloc(2, sizeof(**args->net));
    }
    args->net[0][0] = values.bounds[0];
    for (int i = 1; i < args->count - 1; ++i) {
        args->net[i][0] = args->net[i - 1][0] + args->h;
    }
    args->net[args->count - 1][0] = values.bounds[1];
    edge_problem(&equation, &values, args);
    printf("Решение краевой задачи методом прогонки 2 порядка точности\nListPlot[{");
    for (int i = 0; i < args->count - 1; ++i) {
        printf("{%.10g, %.10g}, ", args->net[i][0], args->net[i][1]);
    }
    printf("{%.10g, %.10g}}]\n", args->net[args->count - 1][0], args->net[args->count - 1][1]);
    for (int i = 0; i < args->count; ++i) {
        free(args->net[i]);
    }
    free(args->net);
    free(args);
    return 0;
}