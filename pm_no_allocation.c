// 第二部分的code，尚未平行

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define pi acos(-1)

// Choose the calculation type
#define NGP 1
#define CLC 2
#define TSC 3
#define type NGP

const double L = 1.0; // 1-D computational domain size
const int N = 5;      // number of equally spaced sampling points
const int np = 10;    // # of particles
const int G = 1;      //gravitational const

double t = 0.0;
double dt = 0.01; // Time scale

// grid length
double dx = L / N;
// Define functions

void linspace(float xi, float xf, int n, double func[]);
void csv_converter(char name[], double arr[np][3], double m[np]);
void set_Narray0_3d(double arr[N][N][N]);
void NGP_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void CLC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void TSC_deposition(double p[np][3], double pm[np], double *g[3], double g_m[N][N][N], int g_sp[np][3], double g_w[np][type][type][type]);
void force_interpolation(double p[np][3], double pm[np], double *g[3], double g_f[N][N][N][3], int g_sp[np][3], double g_w[np][type][type][type], double a[np][3]);

int main(void)
{
    // Store particle position
    double pos[np][3];
    double mass[np];

    // Input the position of particles from csv file
    char filename[] = "test.csv";
    csv_converter(filename, pos, mass);

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

    // initialize some storage array
    double grid_mass[N][N][N];                // Store grid mass value
    set_Narray0_3d(grid_mass);                // Set all # to 0
    int grid_start_pos[np][3];                // The smallest grid id
    double grid_weight[np][type][type][type]; // It's a density function. NGP use 1. CLC use 8. TSC use 27. The order is the same with three for loop from x to z

    double *grid[3] = {x, y, z}; //  Use to load the grid position

    // Deposite particle
    if (type == NGP)
    {
        NGP_deposition(pos, mass, grid, grid_mass, grid_start_pos, grid_weight);
    }
    if (type == CLC)
    {
        CLC_deposition(pos, mass, grid, grid_mass, grid_start_pos, grid_weight);
    }
    if (type == TSC)
    {
        TSC_deposition(pos, mass, grid, grid_mass, grid_start_pos, grid_weight);
    }
    // // check grid mass
    // double count = 0.0;
    // for (int i = 0; i < N; i++)
    // {
    //     // int p_pos[3]; // To store the grid number of the current particle
    //     for (int j = 0; j < N; j++)
    //     {
    //         for (int k = 0; k < N; k++)
    //         {
    //             printf("%1.1f  ", grid_mass[i][j][k]);
    //             count += grid_mass[i][j][k];
    //         }
    //     }
    // }
    // printf("\nTotal mass is %f", count);

    // Calculate the gravitational force with Poisson solver
    //~~~~~~~~~~~~~
    //~~~~~~~~~~~~~
    //~~~~~~~~~~~~~
    //~~~~~~~~~~~~~
    //~~~~~~~~~~~~~
    //~~~~~~~~~~~~~

    // Let the result force is
    double grid_force[N][N][N][3];

    // Accelaration for each particle
    double acc[np][3];

    //  Interpolation
    force_interpolation(pos, mass, grid, grid_force, grid_start_pos, grid_weight, acc);

    // if (type == CLC)
    // {

    // }
    // if (type == TSC)
    // {

    // }

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

void csv_converter(char name[], double arr[np][3], double m[np])
{

    char table[np + 1][5][15];

    FILE *fp = fopen(name, "r"); // fp指向文件头部
    for (int i = 0; i < np + 1; i++)
        for (int j = 0; j < 5; j++)
        {
            fscanf(fp, "%[^,\n]", table[i][j]); // 读到逗号或\n为止
            fgetc(fp);                          // 读取一个字符(逗号或\n)
        }
    fclose(fp);

    // convert it to double and store to the array "pos"

    for (int i = 1; i < np + 1; i++)
    {
        for (int j = 1; j < 5; j++)
        {
            if (j == 1)
            {
                char *c = table[i][j];
                double number = strtod(c, NULL);
                m[i - 1] = number;
                continue;
            }
            char *c = table[i][j];
            double number = strtod(c, NULL);
            arr[i - 1][j - 2] = number;
        }
    }
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
                    count += g_w[i][xn][yn][zn];
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
                }
            }
        }

        // check if pos is correct
        printf("%d.  ", i + 1);
        for (int j = 0; j < 3; j++)
        {
            printf("%f -> %f, ", p[i][j], g[j][p_pos[j]]);
        }
        printf("\n");
    }
}

void force_interpolation(double p[np][3], double pm[np], double *g[3], double g_f[N][N][N][3], int g_sp[np][3], double g_w[np][type][type][type], double a[np][3])
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
