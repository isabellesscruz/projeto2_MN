#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_PONTOS 20
#define PI 3.141592653589793

double calcula_f(double x)
{
    return cos(PI * x);
}
double calcula_df(double x)
{
    return -PI * sin(PI * x);
}
const double X_INI = 0.0;
const double X_FIM = 1.0;

void le_pontos_de_arquivo(const char *nome_arquivo, double *x, double *y, int *num_pontos)
{
    FILE *arquivo = fopen(nome_arquivo, "r");
    if (arquivo == NULL)
    {
        perror("Erro ao abrir o arquivo de entrada");
        exit(EXIT_FAILURE);
    }

    fscanf(arquivo, "%d", num_pontos);

    if (*num_pontos > MAX_PONTOS)
    {
        printf("Erro: O numero de pontos (%d) excede o maximo permitido (%d).\n", *num_pontos, MAX_PONTOS);
        fclose(arquivo);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *num_pontos; ++i)
    {
        fscanf(arquivo, "%lf %lf", &x[i], &y[i]);
    }

    fclose(arquivo);
}

// --- Cálculo do Spline Fixado ---
void calcula_spline_fixado(int num_pontos, double *x, double *y,
                           double *a, double *b, double *c, double *d)
{
    double h[MAX_PONTOS];
    double alpha[MAX_PONTOS];
    double l[MAX_PONTOS];
    double mu[MAX_PONTOS];
    double z[MAX_PONTOS];

    for (int i = 0; i < num_pontos - 1; ++i)
    {
        h[i] = x[i + 1] - x[i];
    }

    alpha[0] = 3.0 * (y[1] - y[0]) / h[0] - 3.0 * calcula_df(x[0]);
    alpha[num_pontos - 1] = 3.0 * calcula_df(x[num_pontos - 1]) - 3.0 * (y[num_pontos - 1] - y[num_pontos - 2]) / h[num_pontos - 2];
    for (int i = 1; i < num_pontos - 1; ++i)
    {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    l[0] = 2.0 * h[0];
    mu[0] = 0.5;
    z[0] = alpha[0] / l[0];
    for (int i = 1; i < num_pontos - 1; ++i)
    {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[num_pontos - 1] = h[num_pontos - 2] * (2.0 - mu[num_pontos - 2]);
    z[num_pontos - 1] = (alpha[num_pontos - 1] - h[num_pontos - 2] * z[num_pontos - 2]) / l[num_pontos - 1];
    c[num_pontos - 1] = z[num_pontos - 1];

    for (int j = num_pontos - 2; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        a[j] = y[j];
    }
}

void mostra_polinomios(FILE *arquivo_saida, int num_pontos, double *x,
                       double *a, double *b, double *c, double *d)
{
    fprintf(arquivo_saida, "\n\n========================================\n");
    fprintf(arquivo_saida, "Polinomios Gerados para cada Intervalo\n");
    fprintf(arquivo_saida, "========================================\n");
    fprintf(arquivo_saida, "Formato: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3\n");

    for (int i = 0; i < num_pontos - 1; ++i)
    {
        fprintf(arquivo_saida, "\n--> Intervalo [%.4f, %.4f]:\n", x[i], x[i + 1]);
        fprintf(arquivo_saida, "    S_%d(x) = %.6f", i, a[i]);

        if (b[i] >= 0)
        {
            fprintf(arquivo_saida, " + %.6f", b[i]);
        }
        else
        {
            fprintf(arquivo_saida, " - %.6f", -b[i]);
        }
        fprintf(arquivo_saida, " * (x - %.4f)", x[i]);

        if (c[i] >= 0)
        {
            fprintf(arquivo_saida, " + %.6f", c[i]);
        }
        else
        {
            fprintf(arquivo_saida, " - %.6f", -c[i]);
        }
        fprintf(arquivo_saida, " * (x - %.4f)^2", x[i]);

        if (d[i] >= 0)
        {
            fprintf(arquivo_saida, " + %.6f", d[i]);
        }
        else
        {
            fprintf(arquivo_saida, " - %.6f", -d[i]);
        }
        fprintf(arquivo_saida, " * (x - %.4f)^3\n", x[i]);
    }
}

double avalia_spline(int num_pontos, double *x, double *a, double *b, double *c, double *d, double xp)
{
    int i;
    for (i = num_pontos - 2; i >= 0; --i)
    {
        if (xp >= x[i])
        {
            break;
        }
    }
    double dx = xp - x[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

double integra_spline(int num_pontos, double *x, double *a, double *b, double *c, double *d)
{
    double integral_total = 0.0;
    for (int i = 0; i < num_pontos - 1; ++i)
    {
        double dx = x[i + 1] - x[i];
        integral_total += a[i] * dx + b[i] * dx * dx / 2.0 + c[i] * dx * dx * dx / 3.0 + d[i] * dx * dx * dx * dx / 4.0;
    }
    return integral_total;
}

void deriva_spline(int num_pontos, double *x, double *a, double *b, double *c, double *d,
                   double xp, double *primeira_deriv, double *segunda_deriv)
{
    int i;
    for (i = num_pontos - 2; i >= 0; --i)
    {
        if (xp >= x[i])
        {
            break;
        }
    }
    double dx = xp - x[i];
    *primeira_deriv = b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;
    *segunda_deriv = 2.0 * c[i] + 6.0 * d[i] * dx;
}

int main()
{
    double pontos_x[MAX_PONTOS];
    double pontos_y[MAX_PONTOS];
    double coef_a[MAX_PONTOS];
    double coef_b[MAX_PONTOS];
    double coef_c[MAX_PONTOS];
    double coef_d[MAX_PONTOS];
    int num_pontos;

    FILE *arquivo_saida = fopen("saida.txt", "w");
    if (arquivo_saida == NULL)
    {
        perror("Erro ao criar o arquivo de saida");
        return EXIT_FAILURE;
    }

    le_pontos_de_arquivo("entrada.txt", pontos_x, pontos_y, &num_pontos);

    fprintf(arquivo_saida, "==================================\n");
    fprintf(arquivo_saida, "Spline Fixado - Dados do Arquivo\n");
    fprintf(arquivo_saida, "==================================\n");
    fprintf(arquivo_saida, "Lendo %d pontos do arquivo 'entrada.txt'.\n", num_pontos);

    calcula_spline_fixado(num_pontos, pontos_x, pontos_y,
                          coef_a, coef_b, coef_c, coef_d);

    mostra_polinomios(arquivo_saida, num_pontos, pontos_x, coef_a, coef_b, coef_c, coef_d);

    fprintf(arquivo_saida, "\n\n========================================\n");
    fprintf(arquivo_saida, "Analise do Spline\n");
    fprintf(arquivo_saida, "========================================\n");

    double integral_calc_spline = integra_spline(num_pontos, pontos_x, coef_a, coef_b, coef_c, coef_d);
    double integral_valor_real = sin(PI * X_FIM) / PI - sin(PI * X_INI) / PI;

    fprintf(arquivo_saida, "\nIntegral do Spline: %.10f\n", integral_calc_spline);
    fprintf(arquivo_saida, "Integral Real:      %.10f\n", integral_valor_real);
    fprintf(arquivo_saida, "Erro na Integral:   %.10e\n", fabs(integral_valor_real - integral_calc_spline));

    double primeira_derivada;
    double segunda_derivada;
    deriva_spline(num_pontos, pontos_x, coef_a, coef_b, coef_c, coef_d,
                  0.5, &primeira_derivada, &segunda_derivada);

    double primeira_deriv_real = -PI * sin(PI * 0.5);
    double segunda_deriv_real = -PI * PI * cos(PI * 0.5);

    fprintf(arquivo_saida, "\nPrimeira Derivada (x=0.5): %.10f (Real: %.10f)\n", primeira_derivada, primeira_deriv_real);
    fprintf(arquivo_saida, "Erro Primeira Derivada:    %.10e\n", fabs(primeira_derivada - primeira_deriv_real));
    fprintf(arquivo_saida, "\nSegunda Derivada (x=0.5):  %.10f (Real: %.10f)\n", segunda_derivada, segunda_deriv_real);
    fprintf(arquivo_saida, "Erro Segunda Derivada:     %.10e\n", fabs(segunda_derivada - segunda_deriv_real));

    fclose(arquivo_saida);

    printf("Processamento concluido. Resultados salvos em 'saida.txt'.\n");

    return 0;
}
