#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#define Lx 1.0 //Length of square domain
#define Ly 1.0 //Width of square domain
#define Wall_Temp 1.0 //Top wall Temperature
#define Nx 128 // 4 for question e
#define Ny 128 // 4 for question e
#define a 0.8 // Relaxation Coefficient for SIP Preconditioned Conjugate Gradient Method
#define pi 3.141592654
#define Max_Iter 500
#define Tol 0.000001

using namespace std;

double dx = Lx / Nx, dy = Ly / Ny;
vector<double> T, d, r_current, r_new, b, A_w, A_s, A_p, A_n, A_e, L_w, L_s, L_p, U_n, U_e;
double res, res_0;

void initialize() 
{
    T.resize(Nx * Ny, 0.0);
    d.resize(Nx * Ny, 0.0);
    r_current.resize(Nx * Ny, 0.0);
    r_new.resize(Nx * Ny, 0.0);
    b.resize(Nx * Ny, 0.0);
    A_w.resize(Nx * Ny, 0.0);
    A_s.resize(Nx * Ny, 0.0);
    A_p.resize(Nx * Ny, 0.0);
    A_n.resize(Nx * Ny, 0.0);
    A_e.resize(Nx * Ny, 0.0);
    L_w.resize(Nx * Ny, 0.0);
    L_s.resize(Nx * Ny, 0.0);
    L_p.resize(Nx * Ny, 0.0);
    U_n.resize(Nx * Ny, 0.0);
    U_e.resize(Nx * Ny, 0.0);
}

void compute_matrix() 
{
    int i;
    for (i = 0; i < Nx * Ny; i ++)
    {
        if (i == 0)
        {
            b[i] = 0.0;
            A_w[i] = 0.0;
            A_s[i] = 0.0;
            A_p[i] = 3.0 * (dy / dx + dx / dy);
            A_n[i] = - dx / dy;
            A_e[i] = - dy / dx;
        }
        else if (i > 0 && i < Ny - 1)
        {
            b[i] = 0.0;
            A_w[i] = 0.0;
            A_s[i] = - dx / dy;
            A_p[i] = 3.0 * dy / dx + 2.0 * dx / dy;
            A_n[i] = - dx / dy;
            A_e[i] = - dy / dx;
        }
        else if (i == Ny - 1)
        {
            b[i] = 2.0 * Wall_Temp  * dx / dy;
            A_w[i] = 0.0;
            A_s[i] = - dx / dy;
            A_p[i] = 3.0 * (dy / dx + dx / dy);
            A_n[i] = 0.0;
            A_e[i] = - dy / dx;
        }
        else if ((i + 1) % Ny == 0.0 && i != Ny - 1 && i != Nx * Ny - 1)
        {
            b[i] = 2.0 * Wall_Temp * dx / dy;
            A_w[i] = - dy / dx;
            A_s[i] = - dx / dy;
            A_p[i] = 2.0 * dy / dx + 3.0 * dx / dy;
            A_n[i] = 0.0;
            A_e[i] = - dy / dx;
        }
        else if (i == Nx * Ny - 1)
        {
            b[i] = 2.0 * Wall_Temp * dx / dy;
            A_w[i] = - dy / dx;
            A_s[i] = - dx / dy;
            A_p[i] = 3.0 * (dy / dx + dx / dy);
            A_n[i] = 0.0;
            A_e[i] = 0.0;
        }
        else if (i > (Nx - 1) * Ny && i < Nx * Ny - 1)
        {
            b[i] = 0.0;
            A_w[i] = - dy / dx;
            A_s[i] = - dx / dy;
            A_p[i] = 3.0 * dy / dx + 2.0 * dx / dy;
            A_n[i] = - dx / dy;
            A_e[i] = 0.0;
        }
        else if (i == (Nx -1) * Ny)
        {
            b[i] = 0.0;
            A_w[i] = - dy / dx;
            A_s[i] = 0.0;
            A_p[i] = 3.0 * (dy / dx + dx / dy);
            A_n[i] = - dx / dy;
            A_e[i] = 0.0;
        }
        else if (i % Ny == 0.0 && i != 0 && i != (Nx - 1) * Ny)
        {
            b[i] = 0.0;
            A_w[i] = - dy / dx;
            A_s[i] = 0.0;
            A_p[i] = 2.0 * dy / dx + 3.0 * dx / dy;
            A_n[i] = - dx / dy;
            A_e[i] = - dy / dx;
        }
        else
        {
            b[i] = 0.0;
            A_w[i] = - dy / dx;
            A_s[i] = - dx / dy;
            A_p[i] = 2.0 * (dy / dx + dx / dy);
            A_n[i] = - dx / dy;
            A_e[i] = - dy / dx;
        }
    }
}

void compute_residual(vector<double>& r)
{
    int i;
    for (i = 0; i < Nx * Ny; i ++)
        {
            if (i == 0)
            {
                r[i] = b[i]- A_p[i] * T[i] - A_n[i] * T[i + 1] - A_e[i] * T[i + Ny];
            }
            else if (i > 0 && i < Ny - 1)
            {
                r[i] = b[i] - A_s[i] * T[i - 1] - A_p[i] * T[i] - A_n[i] * T[i + 1] - A_e[i] * T[i + Ny];
            }
            else if (i == Ny - 1)
            {
                r[i] = b[i]- A_s[i] * T[i - 1] - A_p[i] * T[i] - A_e[i] * T[i + Ny];
            }
            else if ((i + 1) % Ny == 0.0 && i != Ny - 1 && i != Nx * Ny - 1)
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_s[i] * T[i - 1] - A_p[i] * T[i] - A_e[i] * T[i + Ny];
            }
            else if (i == Nx * Ny - 1)
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_s[i] * T[i - 1] - A_p[i] * T[i];
            }
            else if (i > (Nx - 1) * Ny && i < Nx * Ny - 1)
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_s[i] * T[i - 1] - A_p[i] * T[i] - A_n[i] * T[i + 1];
            }
            else if (i == (Nx -1) * Ny)
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_p[i] * T[i] - A_n[i] * T[i + 1];
            }
            else if (i % Ny == 0.0 && i != 0 && i != (Nx - 1) * Ny)
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_p[i] * T[i] - A_n[i] * T[i + 1] - A_e[i] * T[i + Ny];
            }
            else
            {
                r[i] = b[i] - A_w[i] * T[i - Ny] - A_s[i] * T[i - 1] - A_p[i] * T[i] - A_n[i] * T[i + 1] - A_e[i] * T[i + Ny];
            }
        }
}

void compute_A_d(vector<double>& temp)
{
    int i;
    for (i = 0; i < Nx * Ny; i ++)
    {
        if (i == 0)
        {
            temp[i] = A_p[i] * d[i] + A_n[i] * d[i + 1] + A_e[i] * d[i + Ny];
        }
        else if (i > 0 && i < Ny - 1)
        {
            temp[i] = A_s[i] * d[i - 1] + A_p[i] * d[i] + A_n[i] * d[i + 1] + A_e[i] * d[i + Ny];
        }
        else if (i == Ny - 1)
        {
            temp[i] = A_s[i] * d[i - 1] + A_p[i] * d[i] + A_e[i] * d[i + Ny];
        }
        else if ((i + 1) % Ny == 0.0 && i != Ny - 1 && i != Nx * Ny - 1)
        {
            temp[i] = A_w[i] * d[i - Ny] + A_s[i] * d[i - 1] + A_p[i] * d[i] + A_e[i] * d[i + Ny];
        }
        else if (i == Nx * Ny - 1)
        {
            temp[i] = A_w[i] * d[i - Ny] + A_s[i] * d[i - 1] + A_p[i] * d[i];
        }
        else if (i > (Nx - 1) * Ny && i < Nx * Ny - 1)
        {
            temp[i] = A_w[i] * d[i - Ny] + A_s[i] * d[i - 1] + A_p[i] * d[i] + A_n[i] * d[i + 1];
        }
        else if (i == (Nx -1) * Ny)
        {
            temp[i] = A_w[i] * d[i - Ny] + A_p[i] * d[i] + A_n[i] * d[i + 1];
        }
        else if (i % Ny == 0.0 && i != 0 && i != (Nx - 1) * Ny)
        {
            temp[i] = A_w[i] * d[i - Ny] + A_p[i] * d[i] + A_n[i] * d[i + 1] + A_e[i] * d[i + Ny];
        }
        else
        {
            temp[i] = A_w[i] * d[i - Ny] + A_s[i] * d[i - 1] + A_p[i] * d[i] + A_n[i] * d[i + 1] + A_e[i] * d[i + Ny];
        }
    }
}

void compute_L_and_U_ILU()
{
    ofstream fout;
    int i;
    for (i = 0; i < Nx * Ny; i++)
    {
        if (i < Ny)
        {
            L_w[i] = 0.0;
        }
        else
        {
            L_w[i] = A_w[i];
        }
        if (i == 0)
        {
            L_s[i] = 0.0;
            L_p[i] = A_p[i];
        }
        else if (i < Ny)
        {
            L_s[i] = A_s[i];
            L_p[i] = A_p[i] - L_s[i] * U_n[i - 1];
        }
        else
        {
            L_s[i] = A_s[i];
            L_p[i] = A_p[i] - L_s[i] * U_n[i - 1] - L_w[i] * U_e[i - Ny];
        }
        if (i < Nx * Ny - 1)
        {
            U_n[i] = A_n[i] / L_p[i];
        }
        else
        {
            U_n[i] = 0.0;
        }
        if (i < (Nx - 1) * Ny)
        {
            U_e[i] = A_e[i] / L_p[i];
        }
        else
        {
            U_e[i] = 0.0;
        }
    }
    // To print L and U matrices for Ni = Nj = 4
    fout.open("L_U_ILU.dat");
    for (i = 0; i < Nx * Ny; i ++)
    {
        fout << i + 1 << "\t\t\t" << L_w[i] << "\t\t\t" << L_s[i] << "\t\t\t" << L_p[i] << "\t\t\t" << U_n[i] << "\t\t\t" << U_e[i] << endl;
    }
    fout.close();
}

void compute_L_and_U_SIP()
{
    ofstream fout;
    int i;
    for (i = 0; i < Nx * Ny; i++)
    {
        if (i < Ny)
        {
            L_w[i] = 0.0;
        }
        else
        {
            L_w[i] = A_w[i] / (1 +  a * U_n[i - Ny]);
        }
        if (i == 0)
        {
            L_s[i] = 0.0;
            L_p[i] = A_p[i];
        }
        else if (i < Ny)
        {
            L_s[i] = A_s[i] / (1 + a * U_e[i - 1]);
            L_p[i] = A_p[i] + L_s[i] * (a * U_e[i - 1] - U_n[i - 1]);
        }
        else
        {
            L_s[i] = A_s[i] / (1 + a * U_e[i - 1]);
            L_p[i] = A_p[i] + L_w[i] * (a * U_n[i - Ny] - U_e[i - Ny]) + L_s[i] * (a * U_e[i - 1] - U_n[i - 1]);
        }
        if (i < Ny)
        {
            U_n[i] = A_n[i] / L_p[i];
        }
        else if (i < Nx * Ny - 1)
        {
            U_n[i] = (A_n[i] - a * L_w[i] * U_n[i - Ny]) / L_p[i];
        }
        else
        {
            U_n[i] = 0.0;
        }
        if (i == 0)
        {
            U_e[i] = A_e[i] / L_p[i];
        }
        else if (i < (Nx - 1) * Ny)
        {
            U_e[i] = (A_e[i] - a * L_s[i] * U_e[i - 1]) / L_p[i];
        }
        else
        {
            U_e[i] = 0.0;
        }
    }
    // To print L and U matrices for Ni = Nj = 4
    fout.open("L_U_SIP.dat");
    for (i = 0; i < Nx * Ny; i ++)
    {
        fout << i + 1 << "\t\t\t" << L_w[i] << "\t\t\t" << L_s[i] << "\t\t\t" << L_p[i] << "\t\t\t" << U_n[i] << "\t\t\t" << U_e[i] << endl;
    }
    fout.close();
}

void lu_preconditioner(vector<double>& delta, vector<double>& r)
{
    vector<double> R(Nx * Ny, 0.0);
    int i;
    for (i = 0; i < Nx * Ny; i ++)
    {
        if (i == 0)
        {
            R[i] = r[i] / L_p[i];
        }
        else if (i < Ny)
        {
            R[i] = (r[i] - L_s[i] * R[i - 1]) / L_p[i];
        }
        else
        {
            R[i] = (r[i] - L_s[i] * R[i - 1] - L_w[i] * R[i - Ny]) / L_p[i];
        }
    }
    delta[Nx * Ny - 1] = R[Nx * Ny - 1];
    for (i = Nx * Ny - 2; i >= 0; i --)
    {
        if (i >= (Nx - 1) * Ny)
        {
            delta[i] = R[i] - U_n[i] * delta[i + 1];
        }
        else
        {
            delta[i] = R[i] - U_n[i] * delta[i + 1] - U_e[i] * delta[i + Ny];
        }
    }
}

void analytical_solution()
{
    int i, j, n = 1;
    while (n <= 100)
    {
        for (i = 0; i < Nx; i ++)
        {
            for (j = 0; j < Ny; j ++)
            {
                T[i * Nx + j] += 2.0 / pi * (pow(-1, n + 1) + 1) / n * sin(n * pi * (i + 0.5) * dx) * sinh(n * pi * (j + 0.5) * dy) / sinh(n * pi);
            }
        }
        n++;
    }
    
}
void conjugate_gradient(int k)
{
    int i;
    double alpha, alpha_num = 0.0, alpha_den = 0.0, beta, beta_num = 0.0, beta_den = 0.0;
    vector<double> temp(Nx * Ny, 0.0);
    ofstream fout;
    
    if (k == 1)
    {
        res_0 = 0.0;
        compute_residual(r_current);
        for (i = 0; i < Nx * Ny; i ++)
        {
            d[i] = r_current[i];
        }
        for (i = 0; i < Nx * Ny; i ++)
        {
            res_0 += pow(r_current[i], 2);
        }
        res_0 = sqrt(res_0);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_current[i] = r_new[i];
        }
    }
    res = 0.0;
    // Compute alpha
    compute_A_d(temp);
    for (i = 0; i < Nx * Ny; i ++)
    {
        alpha_num += r_current[i] * r_current[i];
        alpha_den += d[i] * temp[i];
    }
    alpha = alpha_num / alpha_den;
    // Compute T and residual
    for (i = 0; i < Nx * Ny; i ++)
    {
        T[i] = T[i] + alpha * d[i];
    }
    if (k % 50 == 0.0)
    {
        compute_residual(r_new);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_new[i] = r_current[i] - alpha * temp[i];
        }
    }
    // Compute beta
    for (i = 0; i < Nx * Ny; i ++)
    {
        beta_num += r_new[i] * r_new[i];
        beta_den += r_current[i] * r_current[i];
    }
    beta = beta_num / beta_den;
    for (i = 0; i < Nx * Ny; i ++)
    {
        d[i] = r_new[i] + beta * d[i];
    }
    for (i = 0; i < Nx * Ny; i ++)
    {
        res += pow(r_new[i], 2);
    }
    res = sqrt(res);
    fout.open("iterations_vs_residual_2_norm_conjugate_gradient.dat", ios::app);
    fout << k << "\t\t" << res << endl;
    fout.close();
}

void jacobi_conjugate_gradient(int k)
{
    int i;
    double alpha, alpha_num = 0.0, alpha_den = 0.0, beta, beta_num = 0.0, beta_den = 0.0;
    vector<double> temp(Nx * Ny, 0.0);
    ofstream fout;
    
    if (k == 1)
    {
        res_0 = 0.0;
        res_0 = 0.0;
        compute_residual(r_current);
        for (i = 0; i < Nx * Ny; i ++)
        {
            d[i] = r_current[i] / A_p[i];
        }
        for (i = 0; i < Nx * Ny; i ++)
        {
            res_0 += pow(r_current[i], 2);
        }
        res_0 = sqrt(res_0);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_current[i] = r_new[i];
        }
    }
    res = 0.0;
    // Compute alpha
    compute_A_d(temp);
    for (i = 0; i < Nx * Ny; i ++)
    {
        alpha_num += r_current[i] * r_current[i] / A_p[i];
        alpha_den += d[i] * temp[i];
    }
    alpha = alpha_num / alpha_den;
    // Compute T and residual
    for (i = 0; i < Nx * Ny; i ++)
    {
        T[i] = T[i] + alpha * d[i];
    }
    if (k % 50 == 0.0)
    {
        compute_residual(r_new);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_new[i] = r_current[i] - alpha * temp[i];
        }
    }
    // Compute beta
    for (i = 0; i < Nx * Ny; i ++)
    {
        beta_num += r_new[i] * r_new[i] / A_p[i];
        beta_den += r_current[i] * r_current[i] / A_p[i];
    }
    beta = beta_num / beta_den;
    for (i = 0; i < Nx * Ny; i ++)
    {
        d[i] = r_new[i] / A_p[i] + beta * d[i];
    }
    for (i = 0; i < Nx * Ny; i ++)
    {
        res += pow(r_new[i], 2);
    }
    res = sqrt(res);
    fout.open("iterations_vs_residual_2_norm_jacobi_conjugate_gradient.dat", ios::app);
    fout << k << "\t\t" << res << endl;
    fout.close();
}

void ilu_conjugate_gradient(int k)
{
    int i;
    double alpha, alpha_num = 0.0, alpha_den = 0.0, beta, beta_num = 0.0, beta_den = 0.0;
    vector<double> temp(Nx * Ny, 0.0), delta_1(Nx * Ny, 0.0), delta_2(Nx * Ny, 0.0);
    ofstream fout;
    
    if (k == 1)
    {
        res_0 = 0.0;
        compute_residual(r_current);
        for (i = 0; i < Nx * Ny; i ++)
        {
            res_0 += pow(r_current[i], 2);
        }
        res_0 = sqrt(res_0);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_current[i] = r_new[i];
        }
    }
    res = 0.0;
    lu_preconditioner(delta_1, r_current);
    if (k == 1)
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            d[i] = delta_1[i];
        }
    }
    // Compute alpha
    compute_A_d(temp);
    for (i = 0; i < Nx * Ny; i ++)
    {
        alpha_num += r_current[i] * delta_1[i];
        alpha_den += d[i] * temp[i];
    }
    alpha = alpha_num / alpha_den;
    // Compute T and residual
    for (i = 0; i < Nx * Ny; i ++)
    {
        T[i] = T[i] + alpha * d[i];
    }
    if (k % 8 == 0.0)
    {
        compute_residual(r_new);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_new[i] = r_current[i] - alpha * temp[i];
        }
    }
    lu_preconditioner(delta_2, r_new);
    // Compute beta
    for (i = 0; i < Nx * Ny; i ++)
    {
        beta_num += r_new[i] * delta_2[i];
        beta_den += r_current[i] * delta_1[i];
    }
    beta = beta_num / beta_den;
    for (i = 0; i < Nx * Ny; i ++)
    {
        d[i] = delta_2[i] + beta * d[i];
    }
    for (i = 0; i < Nx * Ny; i ++)
    {
        res += pow(r_new[i], 2);
    }
    res = sqrt(res);
    fout.open("iterations_vs_residual_2_norm_ilu_conjugate_gradient.dat", ios::app);
    fout << k << "\t\t" << res << endl;
    fout.close();
}

void sip_conjugate_gradient(int k)
{
    int i;
    double alpha, alpha_num = 0.0, alpha_den = 0.0, beta, beta_num = 0.0, beta_den = 0.0;
    vector<double> temp(Nx * Ny, 0.0), delta_1(Nx * Ny, 0.0), delta_2(Nx * Ny, 0.0);
    ofstream fout;
    
    if (k == 1)
    {
        res_0 = 0.0;
        compute_residual(r_current);
        for (i = 0; i < Nx * Ny; i ++)
        {
            res_0 += pow(r_current[i], 2);
        }
        res_0 = sqrt(res_0);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_current[i] = r_new[i];
        }
    }
    res = 0.0;
    lu_preconditioner(delta_1, r_current);
    if (k == 1)
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            d[i] = delta_1[i];
        }
    }
    // Compute alpha
    compute_A_d(temp);
    for (i = 0; i < Nx * Ny; i ++)
    {
        alpha_num += r_current[i] * delta_1[i];
        alpha_den += d[i] * temp[i];
    }
    alpha = alpha_num / alpha_den;
    // Compute T and residual
    for (i = 0; i < Nx * Ny; i ++)
    {
        T[i] = T[i] + alpha * d[i];
    }
    if (k % 8 == 0.0)
    {
        compute_residual(r_new);
    }
    else
    {
        for (i = 0; i < Nx * Ny; i ++)
        {
            r_new[i] = r_current[i] - alpha * temp[i];
        }
    }
    lu_preconditioner(delta_2, r_new);
    // Compute beta
    for (i = 0; i < Nx * Ny; i ++)
    {
        beta_num += r_new[i] * delta_2[i];
        beta_den += r_current[i] * delta_1[i];
    }
    beta = beta_num / beta_den;
    for (i = 0; i < Nx * Ny; i ++)
    {
        d[i] = delta_2[i] + beta * d[i];
    }
    for (i = 0; i < Nx * Ny; i ++)
    {
        res += pow(r_new[i], 2);
    }
    res = sqrt(res);
    fout.open("iterations_vs_residual_2_norm_sip_conjugate_gradient.dat", ios::app);
    fout << k << "\t\t" << res << endl;
    fout.close();
}

void compute_center_temperatures()
{
    int i;
    vector<double>Tx_center(Nx + 1, 0.0);
    vector<double>Ty_center(Ny + 1, 0.0);
    ofstream fout;
    Tx_center[0] = 0.0;
    Tx_center[Nx] = 0.0;
    Ty_center[0] = 0.0;
    Ty_center[Ny] = Wall_Temp;
    for (i = 1; i < Nx; i ++)
    {
        Tx_center[i] = (T[(i - 1 / 2) * Ny - 1] + T[(i - 1 / 2) * Ny] +  T[(i + 1 / 2) * Ny - 1] + T[(i + 1 / 2) * Ny]) / 4.0;
    }
    for (i = 1; i < Ny; i ++)
    {
        Ty_center[i] = (T[Ny * (Nx / 2 - 1) + i - 1] + T[Ny * (Nx / 2 - 1) + i] + T[Ny * Nx / 2 + i - 1] + T[Ny * Nx / 2 + i]) / 4.0;
    }
    fout.open("Temperature_at_x_0.5.dat");
    for (i = 0; i <= Ny; i ++)
    {
        fout << i * dy << "\t\t" << Ty_center[i] << endl;
    }
    fout.close();
    fout.open("Temperature_at_y_0.5.dat");
    for (i = 0; i <= Nx; i ++)
    {
        fout << i * dx << "\t\t" << Tx_center[i] << endl;
    }
    fout.close();
}
void outputdata ()
{
    ofstream fout;
    int i, j;
    fout.open("Temperature.dat");
    for (i = 0; i < Nx; i ++)
    {
        for (j = 0; j < Ny; j ++)
        {
            fout << (i + 0.5) * dx << "\t\t" << (j + 0.5) * dy << "\t\t" << T[i * Nx + j] << endl;
        }
    }
    fout.close();
}
int main()
{
    int count = 1;
    initialize(); // Function to initialize temperature and other data arrays
    //analytical_solution(); // Compute analytical solution
    compute_matrix(); // Compute matrices A and b
    //compute_L_and_U_ILU(); // To be used for ILU method
    compute_L_and_U_SIP(); // To be used for SIP method
    do
    {
        cout << "Iteration no: " << count << endl;
	// Choose any one of the methods that is to be used by removing the double slash before it
        // conjugate_gradient(count);
        // jacobi_conjugate_gradient(count);
        // ilu_conjugate_gradient(count);
        // sip_conjugate_gradient(count);
        if (res < Tol)
        {
            cout << "Termination condition is satisfied" << endl;
            break;
        }
        count++;
    } while (count <= Max_Iter);
    compute_center_temperatures();
    outputdata();
}