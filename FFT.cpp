#include <fftw3.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Choose the deposition type
#define NGP 1
#define CLC 2
#define TSC 3
#define type TSC

// Choose the integration type
#define KDK 1
#define DKD 2
#define type_int KDK

const int np = 1; // # of particles        MUSt BE SET AT FIRST
const int G = 1;  //gravitational const

const double pi = 3.14159265359;
const double L = 2. * pi;
const int N = 8;

int step = 0; // record # of steps
double t = 0.0;
double dt = 0.01; // Time scale

// determine the runing step or total time. You can choose one of them to set
int total_step = 100;
double tf = 1; //end time

// grid length
double dx = L / N;
// Define functions

void linspace(float xi, float xf, int n, double func[]);
void csv_converter(char name[], double p[np][3], double pv[np][3], double m[np]);
void data_output(double p[np][3], double time);
void set_Narray0_3d(double arr[N][N][N]);
void set_particle_param0(double arr[np][3]);
void NGP_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void CLC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void TSC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void force_interpolation(double p[np][3], double pm[np], double g_f[N][N][N][3], int g_sp[np][3], double g_w[np][type][type][type], double a[np][3]);
void K_move(double a[np][3], double pv[np][3], double n);
void D_move(double p[np][3], double pv[np][3], double n);
void get_acc(double p[np][3], double pm[np], double *g[3], double a[np][3]);
void set_3d_vec0(double arr[N][N][N][3]);
void poisson(double grid_mass[][N][N], double grid_acc[][N][N][3])
    // The following functions are just for testing
    void direct_N(double p[np][3], double pm[np], double a[np][3]);
void direct_N_mesh(double g_m[N][N][N], double *g[3], double g_a[N][N][N][3]);

int main(void)
{
    // Store particle position
    double pos[np][3];
    double velocity[np][3];
    double mass[np];

    // Input the position of particles from csv file
    char filename[] = "test.csv";
    csv_converter(filename, pos, velocity, mass);

    // // check if "pos" and "mass" is OK~
    // for (int i = 0; i < np; i++)
    // {
    //     printf("%d.  %f  ", i, mass[i]);

    //     for (int j = 0; j < 3; j++)
    //     {
    //         printf("%f  ", pos[i][j]);
    //     }
    //     printf("\n");
    // }

    // Create 3D coordinates, does not include both ends (grid point)
    double x[N], y[N], z[N];
    linspace(0.0, L, N, x);
    linspace(0.0, L, N, y);
    linspace(0.0, L, N, z);
    double *grid[3] = {x, y, z}; //  Use to load the grid position

    while (t < tf || step < total_step)
    {
        double acc[np][3];
        set_particle_param0(acc);
        if (type_int == KDK)
        {
            get_acc(pos, mass, grid, acc);
            K_move(acc, velocity, 0.5);
            D_move(pos, velocity, 1);
            set_particle_param0(acc);
            get_acc(pos, mass, grid, acc);
            K_move(acc, velocity, 0.5);
        }
        if (type_int == DKD)
        {
            D_move(pos, velocity, 0.5);
            get_acc(pos, mass, grid, acc);
            K_move(acc, velocity, 1);
            D_move(pos, velocity, 0.5);
        }
    }
    return 0;
}

// Use Dynamic Allocation to form a linspace like Python
void linspace(float xi, float xf, int n, double func[])
{
    float space = (xf - xi) / n;
    // printf("%f\n", space);
    for (int i = 0; i < n; i++)
    {
        func[i] = xi + space / 2 + i * space;
    }
}

void csv_converter(char name[], double p[np][3], double pv[np][3], double m[np])
{

    char table[np + 1][5][15];

    FILE *fp = fopen(name, "r"); // fp指向文件头部
    for (int i = 0; i < np + 1; i++)
        for (int j = 0; j < 8; j++)
        {
            fscanf(fp, "%[^,\n]", table[i][j]); // 读到逗号或\n为止
            fgetc(fp);                          // 读取一个字符(逗号或\n)
        }
    fclose(fp);

    // convert it to double and store to the array "pos"

    for (int i = 1; i < np + 1; i++)
    {
        for (int j = 1; j < 8; j++)
        {
            if (j == 1)
            {
                char *c = table[i][j];
                double number = strtod(c, NULL);
                m[i - 1] = number;
                // printf("%f ", number);
            }
            else if (j == 2 || j == 3 || j == 4)
            {
                char *c = table[i][j];
                double number = strtod(c, NULL);
                p[i - 1][j - 2] = number;
                // printf("%f ", number);
            }
            else
            {
                char *c = table[i][j];
                double number = strtod(c, NULL);
                pv[i - 1][j - 5] = number;
                // printf("%f ", number);
            }
        }
        printf("\n");
    }
}

void data_output(double p[np][3], double time)
{
    FILE *fp = fopen("data.txt", "a");
    fprintf(fp, "t = %f\n", time);
    for (int i = 0; i < np; i++)
    {
        fprintf(fp, "%d. (%f, %f, %f)\n", i + 1, p[i][0], p[i][1], p[i][2]);
        printf("%d. (%f, %f, %f)\n", i + 1, p[i][0], p[i][1], p[i][2]);
    }
    fprintf(fp, "----------\n"); // print 10 "-"
}

void set_Narray0_3d(double arr[N][N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                arr[i][j][k] = 0.0;
            }
        }
    }
}

void set_particle_param0(double arr[np][3])
{
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            arr[i][j] = 0.0;
        }
    }
}

void NGP_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type])
{
    for (int i = 0; i < np; i++)
    {
        int p_pos[3]; // To store the grid number of the current particle
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < N; k++)
            {
                if (p[i][j] >= g[j][N - 1])
                {
                    p_pos[j] = N - 1;
                    break;
                }
                if (p[i][j] > g[j][k])
                {
                    continue;
                }
                double l1 = fabs(p[i][j] - g[j][k]);
                if (l1 < dx / 2)
                {
                    p_pos[j] = k;
                    break;
                }
                else
                {
                    p_pos[j] = k - 1;
                    break;
                }
            }
        }

        // add mass
        printf("(%d, %d, %d)\n", p_pos[0], p_pos[1], p_pos[2]);
        g_m[p_pos[0]][p_pos[1]][p_pos[2]] += pm[i];

        // record the dist. starting point
        g_sp[i][0] = p_pos[0];
        g_sp[i][1] = p_pos[1];
        g_sp[i][2] = p_pos[2];

        // record the mass distribution
        g_w[np][0][0][0] = 1;

        // // check if pos is correct
        // printf("%d.  ", i + 1);
        // for (int j = 0; j < 3; j++)
        // {
        //     // if (fabs(g[j][p_pos[j]] - p[i][j]) > dx / 2)
        //     // {
        //     //     printf("error!!(%f > %f)  ", dx / 2, fabs(g[j][p_pos[j]] - p[i][j]));
        //     //     continue;
        //     // }
        //     printf("%f -> %f, ", p[i][j], g[j][p_pos[j]]);
        // }
        // printf("\n");
    }
}

void CLC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type])
{
    for (int i = 0; i < np; i++)
    {
        int p_pos[3]; // To store the grid number of the current particle. It's the start point
        double dist[3][2];
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < N; k++)
            {
                if (p[i][j] >= g[j][N - 1])
                {
                    p_pos[j] = N - 2;
                    dist[j][0] = 0;
                    dist[j][1] = 1;
                    break;
                }
                if (p[i][j] > g[j][k])
                {
                    continue;
                }
                if (p[i][j] <= g[j][0])
                {
                    p_pos[j] = 0;
                    dist[j][0] = 1;
                    dist[j][1] = 0;
                    break;
                }
                double l1 = fabs(p[i][j] - g[j][k]); // In fact, it equals g[j][k] - p[i][j]
                p_pos[j] = k - 1;
                dist[j][0] = l1 / dx;
                dist[j][1] = (dx - l1) / dx;
                break;
            }
        }
        // // check starting position
        // printf("(%d, %d, %d)\n", p_pos[0], p_pos[1], p_pos[2]);

        // record the dist. starting point
        g_sp[i][0] = p_pos[0];
        g_sp[i][1] = p_pos[1];
        g_sp[i][2] = p_pos[2];

        // record the mass distribution
        printf("(%d, %d, %d) :\n", p_pos[0], p_pos[1], p_pos[2]);
        double count = 0;
        for (int xn = 0; xn < 2; xn++)
        {
            for (int yn = 0; yn < 2; yn++)
            {
                for (int zn = 0; zn < 2; zn++)
                {
                    g_w[i][xn][yn][zn] = dist[0][xn] * dist[1][yn] * dist[2][zn];
                    printf("(%d, %d, %d, %f)\n", xn, yn, zn, g_w[i][xn][yn][zn]);
                    count += g_w[i][xn][yn][zn];
                }
            }
        }
        printf("count is %f\n", count);

        // add mass
        for (int xn = 0; xn < 2; xn++)
        {
            for (int yn = 0; yn < 2; yn++)
            {
                for (int zn = 0; zn < 2; zn++)
                {
                    g_m[p_pos[0] + xn][p_pos[1] + yn][p_pos[2] + zn] += pm[i] * g_w[i][xn][yn][zn];
                }
            }
        }

        // // check if pos is correct
        // printf("%d.  ", i + 1);
        // for (int j = 0; j < 3; j++)
        // {
        //     printf("%f -> %f, ", p[i][j], g[j][p_pos[j]]);
        // }
        // printf("\n");
    }
}

void TSC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type])
{
    // printf("dx = %f\n", dx);
    for (int i = 0; i < np; i++)
    {
        int p_pos[3]; // To store the grid number of the current particle. It's the start point
        double dist[3][3];
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < N; k++)
            {
                if (p[i][j] >= g[j][N - 2])
                {
                    p_pos[j] = N - 3;
                    double l1 = fabs(p[i][j] - g[j][N - 2]);
                    if (l1 >= dx / 2)
                    {
                        dist[j][0] = 0;
                        dist[j][1] = pow((1.5 * dx - l1) / dx, 2) / 2;
                        dist[j][2] = 1 - pow((1.5 * dx - l1) / dx, 2) / 2;
                        break;
                    }
                    dist[j][0] = pow((0.5 * dx - l1) / dx, 2) / 2;
                    dist[j][1] = 1 - pow((0.5 * dx - l1) / dx, 2) / 2 - pow((0.5 * dx + l1) / dx, 2) / 2;
                    dist[j][2] = pow((0.5 * dx + l1) / dx, 2) / 2;
                    break;
                }
                if (p[i][j] > g[j][k])
                {
                    continue;
                }
                if (p[i][j] <= g[j][1])
                {
                    p_pos[j] = 0;
                    double l1 = fabs(p[i][j] - g[j][1]); // Use 0-th grid point to define
                    if (l1 >= dx / 2)
                    {
                        dist[j][0] = 1 - pow((1.5 * dx - l1) / dx, 2) / 2;
                        dist[j][1] = pow((1.5 * dx - l1) / dx, 2) / 2;
                        dist[j][2] = 0;
                        break;
                    }
                    dist[j][0] = pow((0.5 * dx + l1) / dx, 2) / 2;
                    dist[j][1] = 1 - pow((0.5 * dx + l1) / dx, 2) / 2 - pow((0.5 * dx - l1) / dx, 2) / 2;
                    dist[j][2] = pow((0.5 * dx - l1) / dx, 2) / 2;
                    break;
                }
                double l1 = fabs(p[i][j] - g[j][k]); // In fact, it equals g[j][k] - p[i][j]
                if (l1 <= dx / 2)
                {
                    p_pos[j] = k - 1;
                    dist[j][0] = pow((0.5 * dx + l1) / dx, 2) / 2;
                    dist[j][1] = 1 - pow((0.5 * dx + l1) / dx, 2) / 2 - pow((0.5 * dx - l1) / dx, 2) / 2;
                    dist[j][2] = pow((0.5 * dx - l1) / dx, 2) / 2;
                    break;
                }
                else
                {
                    p_pos[j] = k - 2;
                    dist[j][0] = pow((l1 - 0.5 * dx) / dx, 2) / 2;
                    dist[j][1] = 1 - pow((l1 - 0.5 * dx) / dx, 2) / 2 - pow((1.5 * dx - l1) / dx, 2) / 2;
                    dist[j][2] = pow((1.5 * dx - l1) / dx, 2) / 2;
                    break;
                }
            }
        }
        // // check starting position
        // printf("(%d, %d, %d)\n", p_pos[0], p_pos[1], p_pos[2]);

        // record the dist. starting point
        g_sp[i][0] = p_pos[0];
        g_sp[i][1] = p_pos[1];
        g_sp[i][2] = p_pos[2];

        // record the mass distribution
        // printf("(%d, %d, %d) :\n", p_pos[0], p_pos[1], p_pos[2]);
        double count = 0;
        for (int xn = 0; xn < 3; xn++)
        {
            for (int yn = 0; yn < 3; yn++)
            {
                // printf("%f  ", dist[xn][yn]);
                for (int zn = 0; zn < 3; zn++)
                {
                    g_w[i][xn][yn][zn] = dist[0][xn] * dist[1][yn] * dist[2][zn];
                    // printf("(%d, %d, %d, %f)\n", xn, yn, zn, g_w[i][xn][yn][zn]);
                    // count += g_w[i][xn][yn][zn];
                }
            }
            // printf("\n");
        }
        // printf("count is %f\n", count);

        // add mass
        for (int xn = 0; xn < 3; xn++)
        {
            for (int yn = 0; yn < 3; yn++)
            {
                for (int zn = 0; zn < 3; zn++)
                {
                    g_m[p_pos[0] + xn][p_pos[1] + yn][p_pos[2] + zn] += pm[i] * g_w[i][xn][yn][zn];
                    count += pm[i] * g_w[i][xn][yn][zn];
                    printf("%f ", pm[i] * g_w[i][xn][yn][zn]);
                }
            }
        }
        printf("\nCount = %f\n", count);

        // check if pos is correct
        // printf("%d.  ", i + 1);
        // for (int j = 0; j < 3; j++)
        // {
        //     printf("%f -> %f, ", p[i][j], g[j][p_pos[j]]);
        // }
        // printf("\n");
    }
}

void force_interpolation(double p[np][3], double pm[np], double g_f[N][N][N][3], int g_sp[np][3], double g_w[np][type][type][type], double a[np][3])
{
    for (int i = 0; i < np; i++)
    {
        for (int xn = 0; xn < type; xn++)
        {
            for (int yn = 0; yn < type; yn++)
            {
                for (int zn = 0; zn < type; zn++)
                {
                    for (int cor = 0; cor < 3; cor++)
                    {
                        a[i][cor] = g_f[g_sp[i][0] + xn][g_sp[i][1] + yn][g_sp[i][2] + zn][cor] * g_w[i][xn][yn][zn] / pm[i];
                    }
                }
            }
        }
    }
}

// orbit integration

void K_move(double a[np][3], double pv[np][3], double n)
{
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            pv[i][j] += a[i][j] * n * dt;
        }
    }
}

void D_move(double p[np][3], double pv[np][3], double n)
{
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            p[i][j] += pv[i][j] * n * dt;
        }
    }
}

void get_acc(double p[np][3], double pm[np], double *g[3], double a[np][3])
{

    // initialize some storage array
    double grid_mass[N][N][N];                // Store grid mass value
    set_Narray0_3d(grid_mass);                // Set all # to 0
    int grid_start_pos[np][3];                // The smallest grid id
    double grid_weight[np][type][type][type]; // It's a density function. NGP use 1. CLC use 8. TSC use 27. The order is the same with three for loop from x to z

    // Let the result force is (That is the result from Poisson solver)
    double grid_acc[N][N][N][3];
    set_3d_vec0(grid_acc);

    // Deposite particle
    if (type == NGP)
    {
        NGP_deposition(p, pm, g, grid_mass, grid_start_pos, grid_weight);
    }
    if (type == CLC)
    {
        CLC_deposition(p, pm, g, grid_mass, grid_start_pos, grid_weight);
    }
    if (type == TSC)
    {
        TSC_deposition(p, pm, g, grid_mass, grid_start_pos, grid_weight);
    }

    // Calculate the gravitational force with Poisson solver
    poisson(grid_mass, grid_acc);
    // Use direct N-body to calculate the gravitational force for comparison

    // direct_N_mesh(grid_mass, g, grid_acc);

    //  force interpolation for all particles 加權(?
    force_interpolation(p, pm, grid_acc, grid_start_pos, grid_weight, a);
}

void set_3d_vec0(double arr[N][N][N][3])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                for (int cor = 0; cor < N; cor++)
                {
                    arr[i][j][k][cor] = 0.0;
                }
            }
        }
    }
}

//
// The following functions are just for testing.
//
// If we need to use the direct N-body to calculate the self-gravity on each partticle
void direct_N(double p[np][3], double pm[np], double a[np][3])
{
    for (int i = 0; i < np; i++)
    {
        for (int j = 0; j < np; j++)
        {
            if (i == j)
            {
                continue;
            }
            for (int cor = 0; cor < 3; cor++)
            {
                a[i][cor] += G * pm[j] * (p[j][cor] - p[i][cor]) / pow(fabs(p[j][cor] - p[i][cor]), 3);
            }
        }
    }
}

// Use direct N-body to calculate the self gravity on mesh
void direct_N_mesh(double g_m[N][N][N], double *g[3], double g_a[N][N][N][3])
{
    for (int x1 = 0; x1 < N; x1++)
    {
        for (int y1 = 0; y1 < N; y1++)
        {
            for (int z1 = 0; z1 < N; z1++)
            {
                for (int x2 = 0; x2 < N; x2++)
                {
                    for (int y2 = 0; y2 < N; y2++)
                    {
                        for (int z2 = 0; z2 < N; z2++)
                        {
                            if (x1 == x2 && y1 == y2 && z1 == z2)
                            {
                                continue;
                            }
                            g_a[x1][y1][z1][0] += G * g_m[x2][y2][z2] * (g[0][x2] - g[0][x1]) / pow(fabs(g[0][x2] - g[0][x1]), 3);
                            g_a[x1][y1][z1][1] += G * g_m[x2][y2][z2] * (g[1][x2] - g[1][x1]) / pow(fabs(g[1][x2] - g[1][x1]), 3);
                            g_a[x1][y1][z1][2] += G * g_m[x2][y2][z2] * (g[2][x2] - g[2][x1]) / pow(fabs(g[2][x2] - g[2][x1]), 3);
                        }
                    }
                }
            }
        }
    }
}

void poisson(double grid_mass[][N][N], double grid_acc[][N][N][3])
{

    int i, j, k;
    double *X, *Y, *Z;
    X = (double *)malloc(N * sizeof(double));
    Y = (double *)malloc(N * sizeof(double));
    Z = (double *)malloc(N * sizeof(double));
    for (i = 0; i < N; i++)
    {
        X[i] = -pi + i * dx;
        for (j = 0; j < N; j++)
        {
            Y[j] = -pi + j * dx;
            for (k = 0; k < N; k++)
            {
                Z[k] = -pi + k * dx;
            }
        }
    }

    fftw_complex *out1, *in2, *out2, *in1;
    fftw_plan p1, p2;

    in1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N * N * N));
    out2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N * N * N));
    out1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N * N * N));
    in2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N * N * N));

    p1 = fftw_plan_dft_3d(N, N, N, in1, out1, FFTW_FORWARD, FFTW_MEASURE);
    p2 = fftw_plan_dft_3d(N, N, N, in2, out2, FFTW_BACKWARD, FFTW_MEASURE);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                in1[i * N * N + j * N + k][0] = grid_mass[i][j][k]; // row major ordering
                in1[i * N * N + j * N + k][1] = 0;
            }
        }
    }

    fftw_execute(p1); // FFT forward

    for (i = 0; i < N; i++)
    { // f = g / ( kx^2 + ky^2 + kz^2 )
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                double ksquared = 0;
                in2[i * N * N + j * N + k][0] = 0;
                in2[i * N * N + j * N + k][1] = 0;
                if (2 * i < N)
                {
                    ksquared = ((double)i * i);
                }
                else
                {
                    ksquared = ((double)(N - i) * (N - i));
                }
                if (2 * j < N)
                {
                    ksquared += ((double)j * j);
                }
                else
                {
                    ksquared += ((double)(N - j) * (N - j));
                }
                if (2 * k < N)
                {
                    ksquared += ((double)k * k);
                }
                else
                {
                    ksquared += ((double)(N - k) * (N - k));
                }
                if (ksquared != 0)
                {
                    in2[i * N * N + j * N + k][0] = out1[i * N * N + j * N + k][0] / ksquared;
                    in2[i * N * N + j * N + k][1] = out1[i * N * N + j * N + k][1] / ksquared;
                }
                else
                {
                    in2[i * N * N + j * N + k][0] = 0;
                    in2[i * N * N + j * N + k][1] = 0;
                }
            }
        }
    }

    fftw_execute(p2); //FFT backward

    // checking the results computed

    double phi[N][N][N];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                phi[i][j][k] = out2[i * N * N + j * N + k][0] / (N * N * N);
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                if (i == 0)
                {
                    grid_acc[i][j][k][0] = -(phi[1][j][k] - phi[N - 1][j][k]) / (2. * dx);
                }
                else if (i == (N - 1))
                {
                    grid_acc[i][j][k][0] = -(phi[0][j][k] - phi[N - 2][j][k]) / (2. * dx);
                }
                else
                {
                    grid_acc[i][j][k][0] = -(phi[i + 1][j][k] - phi[i - 1][j][k]) / (2. * dx);
                }

                if (j == 0)
                {
                    grid_acc[i][j][k][1] = -(phi[i][1][k] - phi[i][N - 1][k]) / (2. * dx);
                }
                else if (j == (N - 1))
                {
                    grid_acc[i][j][k][1] = -(phi[i][0][k] - phi[i][N - 2][k]) / (2. * dx);
                }
                else
                {
                    grid_acc[i][j][k][1] = -(phi[i][j + 1][k] - phi[i][j - 1][k]) / (2. * dx);
                }

                if (k == 0)
                {
                    grid_acc[i][j][k][2] = -(phi[i][j][1] - phi[i][j][N - 1]) / (2. * dx);
                }
                else if (i == (N - 1))
                {
                    grid_acc[i][j][k][2] = -(phi[i][j][0] - phi[i][j][N - 2]) / (2. * dx);
                }
                else
                {
                    grid_acc[i][j][k][2] = -(phi[i][j][k + 1] - phi[i][j][k - 1]) / (2. * dx);
                }
            }
        }
    }
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(out1);
    fftw_free(out2);
    fftw_free(in1);
    fftw_free(in2);
}
