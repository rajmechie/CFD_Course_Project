#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#define Lx 1.0 //Length of square cavity
#define Ly 1.0 //Width of square cavity
#define Lid_vel 1.0 //Top lid velocity
#define Nx 128
#define Ny 128
#define Dt 0.003
#define Max_Iter 10000
#define Max_Timestep 10000
#define Tol 0.00001

// Ghost cells are used on both sides to apply the boundary conditions. As the cell faces touch the cavity walls,
// the boundary conditions are modified such that the wall velocities, pressures and stream function values are satisfied.

double dx = Lx / Nx, dy = Ly / Ny;
double u[Nx + 2][Ny + 2], v[Nx + 2][Ny + 2], u_old[Nx + 2][Ny + 2], v_old[Nx + 2][Ny + 2], p[Nx + 2][Ny + 2], vel_magn[Nx + 2][Ny + 2], omega[Nx + 2][Ny + 2], psi[Nx + 2][Ny + 2];
double u_center[Ny + 1], v_center[Nx + 1];

// Cu_n1, Cu_n2 are used to store the convective terms in u-momentum equation for nth and (n-1)th timestep
// Cv_n1, Cv_n2 are used to store the convective terms in v-momentum equation for nth and (n-1)th timestep
// Du, Dv are used to store the diffusion terms in u and v momentum equation for nth timestep

double Cu_n1[Nx + 2][Ny + 2], Cu_n2[Nx + 2][Ny + 2], Cv_n1[Nx + 2][Ny + 2], Cv_n2[Nx + 2][Ny + 2], 
       Du[Nx + 2][Ny + 2], Dv[Nx + 2][Ny + 2];

using namespace std;

void initialize() 
{
    int i, j;
    for (i = 0; i <= Nx + 1; i ++) 
    {
        for (j = 0; j <= Ny + 1; j ++) 
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            u_old[i][j] = 0.0;
            v_old[i][j] = 0.0;
            p[i][j] = 0.0;
            vel_magn[i][j] = 0.0;
            omega[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
    }
}

void apply_boundary_conditions() 
{
    int i, j;
    for (i = 1; i <= Nx; i ++) 
    {
        u[i][0] = 2 * Lid_vel - u[i][1];  
        v[i][0] = - v[i][1];

        u[i][Ny + 1] = -u[i][Ny]; 
        v[i][Ny + 1] = -v[i][Ny];
    }

    for (j = 1; j <= Ny; j ++) 
    {
        u[0][j] = -u[1][j];    
        v[0][j] = -v[1][j]; 

        u[Nx + 1][j] = -u[Nx][j];  
        v[Nx + 1][j] = -v[Nx][j];
    }
}

void thomas_algorithm(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d,
                      vector<double>& f) 
{
    vector<double> c_star(Ny + 2, 0.0);
    vector<double> d_star(Ny + 2, 0.0);
    int i;
    double m;
                                                                                                                                                           
    c_star[1] =  - c[1] / b[1];
    d_star[1] = d[1] / b[1] - (a[1] / b[1]) * f[0];
                                                                                                                                                     
    for (i = 2; i <= Ny; i ++)
    {
        m = 1.0 / (b[i] + a[i] * c_star[i-1]);
        c_star[i] = - c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
    }

    for (i = Ny; i >= 1; i--) 
    {
    f[i] = d_star[i] + c_star[i] * f[i + 1];
    }                                                                                                                                                                   
}


void compute_u_velocities(double Re)
{
    vector<double> a(Ny + 2, 0.0);
    vector<double> b(Ny + 2, 0.0);
    vector<double> c(Ny + 2, 0.0);
    vector<double> d(Ny + 2, 0.0);
    vector<double> f(Ny + 2, 0.0);
    
    double U, error, error_max;
    int i,j, count = 0;
    do
    {
        error_max = 0.0;
        for (i = 1; i <= Nx; i ++)
        {
            a[0] = 1.0;
            a[Ny + 1] = 1.0;
            b[0] = 0.0;
            b[Ny + 1] = 0.0;
            c[0] = 0.0;
            c[Ny + 1] = 0.0;
            d[0] = u[i][0];
            d[Ny + 1] = u[i][Ny + 1];
            f[0] = u[i][0];
            f[Ny + 1] = u[i][Ny + 1];
            for (j = 1; j <= Ny; j ++)
            {
                a[j] = - Dt / (2.0 * Re * pow(dy, 2));
                b[j] = 1 + Dt / (Re * pow(dx, 2)) + Dt / (Re * pow(dy, 2));
                c[j] = - Dt / (2.0 * Re * pow(dy, 2));
                d[j] = u_old[i][j] - (Dt / 2) * (3 * Cu_n1[i][j] - Cu_n2[i][j]) + (Dt / 2) * Du[i][j] + 
                            (Dt / (2.0 * Re * pow(dx, 2))) * u[i+1][j] + (Dt / (2.0 * Re * pow(dx, 2))) * u[i-1][j];
            }
            thomas_algorithm(a, b, c, d, f);
            for (j = 1; j <= Ny; j ++)
            {
                U = u[i][j];
                u[i][j] = f[j];
                error = fabs(u[i][j] - U);
                if (error >= error_max)
                {
                    error_max = error;
                }
            }
        }
        count ++;
        for (i = 1; i <= Nx; i ++) 
        {
            u[i][0] = 2 * Lid_vel - u[i][1];
            u[i][Ny + 1] = -u[i][Ny];
        }

        for (j = 1; j <= Ny; j ++) 
        {
            u[0][j] = - u[1][j];    
            u[Nx + 1][j] = - u[Nx][j];  
        }
    } while (error_max >= Tol && count <= Max_Iter);

}

void compute_v_velocities(double Re)
{
    vector<double> a(Ny + 2, 0.0);
    vector<double> b(Ny + 2, 0.0);
    vector<double> c(Ny + 2, 0.0);
    vector<double> d(Ny + 2, 0.0);
    vector<double> f(Ny + 2, 0.0);
    
    double V, error, error_max;
    int i,j, count = 0;
    
    do
    {
        error_max =0.0;
        for (i = 1; i <= Nx; i ++)
        {
            a[0] = 1.0;
            a[Ny + 1] = 1.0;
            b[0] = 0.0;
            b[Ny + 1] = 0.0;
            c[0] = 0.0;
            c[Ny + 1] = 0.0;
            d[0] = v[i][0];
            d[Ny + 1] = v[i][Ny + 1];
            f[0] = v[i][0];
            f[Ny + 1] = v[i][Ny + 1];
            for (j = 1; j <= Ny; j ++)
            {
                a[j] = - Dt / (2.0 * Re * pow(dy, 2));
                b[j] = 1 + Dt / (Re * pow(dx, 2)) + Dt / (Re * pow(dy, 2));
                c[j] = - Dt / (2.0 * Re * pow(dy, 2));
                d[j] = v_old[i][j] - (Dt / 2) * (3 * Cv_n1[i][j] - Cv_n2[i][j]) + (Dt / 2) * Dv[i][j] + 
                       (Dt / (2.0 * Re * pow(dx, 2))) * v[i+1][j] + (Dt / (2.0 * Re * pow(dx, 2))) * v[i-1][j];
            }
            thomas_algorithm(a, b, c, d, f);
            for (j = 1; j <= Ny; j ++)
            {
                V = v[i][j];
                v[i][j] = f[j];
                error = fabs(v[i][j] - V);
                if (error >= error_max)
                {
                    error_max = error;
                }
            }
        }
        count ++;
        for (i = 1; i <= Nx; i ++) 
        {
            v[i][0] = - v[i][1];
            v[i][Ny + 1] = - v[i][Ny];
        }

        for (j = 1; j <= Ny; j ++) 
        {
            v[0][j] = -v[1][j]; 
            v[Nx + 1][j] = -v[Nx][j];
        }
    } while (error_max >= Tol && count <= Max_Iter);

}

void solve_momentum_equations(int count, double Re)
{
    double U_west, U_east, V_south, V_north;
    int i,j;
    
    // Computing the convection and diffusion terms
    
    for (j = 1; j <= Ny; j ++)
    {
        for (i = 1; i <= Nx; i ++)
        {
            U_west = (u[i - 1][j] + u[i][j]) / 2;
            U_east = (u[i + 1][j] + u[i][j]) / 2;
            V_south = (v[i][j - 1] + v[i][j]) / 2;
            V_north = (v[i][j + 1] + v[i][j]) / 2;
            if (count != 1)
            {
                Cu_n2[i][j] = Cu_n1[i][j];
                Cv_n2[i][j] = Cv_n1[i][j];
            }
            Cu_n1[i][j] = ((u[i][j] + u[i + 1][j]) * U_east - (u[i - 1][j] + u[i][j]) * U_west) / (2 * dx)
                          + ((u[i][j] + u[i][j + 1]) * V_north - (u[i][j - 1] + u[i][j]) * V_south) / (2 * dy);
            Cv_n1[i][j] = ((v[i][j] + v[i + 1][j]) * U_east - (v[i - 1][j] + v[i][j]) * U_west) / (2 * dx)
                          + ((v[i][j] + v[i][j + 1]) * V_north - (v[i][j - 1] + v[i][j]) * V_south) / (2 * dy);
            if (count == 1)
            {
                Cu_n2[i][j] = Cu_n1[i][j];
                Cv_n2[i][j] = Cv_n1[i][j];
            }
            Du[i][j] = 1 / (Re * pow(dx,2)) * u[i + 1][j] + 1 / (Re * pow(dx,2)) * u[i - 1][j] + 1 / (Re * pow(dy,2)) * u[i][j + 1]
                        + 1 / (Re * pow(dy,2)) * u[i][j - 1] - 2 * (1 / (Re * pow(dx,2)) + 1 / (Re * pow(dy,2))) * u[i][j];
            Dv[i][j] = 1 / (Re * pow(dx,2)) * v[i + 1][j] + 1 / (Re * pow(dx,2)) * v[i - 1][j] + 1 / (Re * pow(dy,2)) * v[i][j + 1]
                        + 1 / (Re * pow(dy,2)) * v[i][j - 1] - 2 * (1 / (Re * pow(dx,2)) + 1 / (Re * pow(dy,2))) * v[i][j];
        }
    }
    
    // Using Line Gauss-Siedel Method
    
    compute_u_velocities(Re);
    compute_v_velocities(Re);
}

void solve_pressure_poisson_equation(double alpha)
{
    vector<double> a(Ny + 2, 0.0);
    vector<double> b(Ny + 2, 0.0);
    vector<double> c(Ny + 2, 0.0);
    vector<double> d(Ny + 2, 0.0);
    vector<double> f(Ny + 2, 0.0);
    
    double U_west, U_east, V_south, V_north, P, error, error_max;
    int i, j, count = 0;
    
    do
    {
        error_max = 0.0;
        for (i = 1; i <= Nx; i ++) 
        {
            p[i][0] = p[i][1];
            p[i][Ny + 1] = p[i][Ny];
        }

        for (j = 1; j <= Ny; j ++) 
        {
            p[0][j] = p[1][j];
            p[Nx + 1][j] = p[Nx][j];
        }
        for (i = 1; i <= Nx; i ++)
        {
            a[0] = 1.0;
            a[Ny + 1] = 1.0;
            b[0] = 0.0;
            b[Ny + 1] = 0.0;
            c[0] = 0.0;
            c[Ny + 1] = 0.0;
            d[0] = p[i][0];
            d[Ny + 1] = p[i][Ny + 1];
            f[0] = p[i][0];
            f[Ny + 1] = p[i][Ny + 1];
            for (j = 1; j <= Ny; j ++)
            {
                U_west = (u[i - 1][j] + u[i][j]) / 2;
                U_east = (u[i + 1][j] + u[i][j]) / 2;
                V_south = (v[i][j - 1] + v[i][j]) / 2;
                V_north = (v[i][j + 1] + v[i][j]) / 2;
                a[j] = - 1.0 / pow(dy, 2);
                b[j] =  2.0 / pow(dx, 2) + 2.0 / pow(dy,2);
                c[j] = - 1.0 / pow(dy, 2);
                d[j] =  - ((U_east - U_west) / dx + (V_north - V_south) / dy) / Dt + p[i + 1][j] / pow(dx, 2) + p[i - 1][j] / pow(dx, 2);
            }
            thomas_algorithm(a, b, c, d, f);
            for (j = 1; j <= Ny; j ++)
            {
                P = p[i][j];
                p[i][j] = P + alpha * (f[j] - P);
                error = fabs(p[i][j] - P);
                if (error >= error_max)
                {
                    error_max = error;
                }
            }
        }
        count ++;
    } while (error_max >= Tol && count <= Max_Iter);
    cout << count << endl;
}

void update_velocity()
{
   int i, j;
   for (i = 1; i <= Nx; i ++)
    {
        for (j = 1; j <= Ny; j ++)
        {
            u[i][j] = u[i][j] - Dt * (p[i + 1][j] - p[i - 1][j]) / (2 * dx);
            v[i][j] = v[i][j] - Dt * (p[i][j + 1] - p[i][j - 1]) / (2 * dy);
        }
    }
}

void compute_velocity_magnitude()
{
    int i, j;
    for (i = 1; i <= Nx; i ++)
        {
            for (j = 1; j <= Ny; j ++)
            {
            vel_magn[i][j] = sqrt(pow(u[i][j], 2) + pow(v[i][j], 2));
            }
        }
}

void compute_vorticity()
{
    int i, j;
    for (i = 1; i <= Nx; i ++)
        {
            for (j = 1; j <= Ny; j ++)
            {
            omega[i][j] = (v[i + 1][j] - v[i - 1][j]) / (2 * dx) - (u[i][j + 1] - u[i][j - 1]) / (2 * dy);
            }
        }
}
void compute_stream_function()
{
    vector<double> a(Ny + 2, 0.0);
    vector<double> b(Ny + 2, 0.0);
    vector<double> c(Ny + 2, 0.0);
    vector<double> d(Ny + 2, 0.0);
    vector<double> f(Ny + 2, 0.0);
    
    double Psi, error, error_max;
    int i, j, count = 0;
    
    do
    {
        error_max = 0.0;
        for (i = 1; i <= Nx; i ++) 
    {
        psi[i][0] = psi[i][1] - dy * Lid_vel;
        psi[i][Ny + 1] = - psi[i][Ny];
    }

    for (j = 1; j <= Ny; j ++) 
    {
        psi[0][j] = - psi[1][j];
        psi[Nx + 1][j] = - psi[Nx][j];
    }
        for (i = 1; i <= Nx; i ++)
        {
            a[0] = 1.0;
            a[Ny + 1] = 1.0;
            b[0] = 0.0;
            b[Ny + 1] = 0.0;
            c[0] = 0.0;
            c[Ny + 1] = 0.0;
            d[0] = psi[i][0];
            d[Ny + 1] = psi[i][Ny + 1];
            f[0] = psi[i][0];
            f[Ny + 1] = psi[i][Ny + 1];
            for (j = 1; j <= Ny; j ++)
            {
                a[j] = - 1.0 / pow(dy, 2);
                b[j] =  2.0 / pow(dx, 2) + 2.0 / pow(dy,2);
                c[j] = - 1.0 / pow(dy, 2);
                d[j] = omega[i][j] + psi[i + 1][j] / pow(dx, 2) + psi[i - 1][j] / pow(dx, 2);
            }
            thomas_algorithm(a, b, c, d, f);
            for (j = 1; j <= Ny; j ++)
            {
                Psi = psi[i][j];
                psi[i][j] = f[j];
                error = fabs(psi[i][j] - Psi);
                if (error >= error_max)
                {
                    error_max = error;
                }
            }
        }
        count ++;
    } while (error_max >= Tol && count <= Max_Iter);
}

void compute_min_streamline()
{
    double psi_min = 0.0, x, y, omega_min;
    int i, j;
    for (i = 1; i <= Nx; i ++)
    {
        for (j = 1; j <= Ny; j ++)
        {
            if (psi[i][j] < psi_min)
            {
                psi_min = psi[i][j];
                x = (i - 0.5) * dx;
                y = Ly - (j - 0.5) * dy;
                omega_min = omega[i][j];
            }
        }
    }
    cout << "Minimum streamline value= " << psi_min << " at x= " << x << ", y= " << y << endl;
    cout << "Vorticity value= " << omega_min << endl;
}

void compute_center_velocities()
{
    int i, j, k;
    u_center[Ny] = 0.0;
    u_center[0] = Lid_vel;
    v_center[0] = 0.0;
    v_center[Nx] = 0.0;
    for (k = 1; k < Ny; k ++)
    {
        u_center[k] = (u[Nx / 2][k] + u[Nx / 2 + 1][k]) / 2;
    }
    for (k = 1; k < Nx; k ++)
    {
        v_center[k] = - (v[k][Ny / 2] + v[k][Ny / 2 + 1]) / 2;
    }
}

void outputdata ()
{
    ofstream fout;
    int i, j;
    fout.open("velocity.dat");
    for (i = 1; i <= Nx; i ++)
    {
        for (j = 1; j <= Ny; j ++)
        {
            fout << (i - 0.5) * dx << "\t\t" << Ly - (j - 0.5) * dy << "\t\t" << vel_magn[i][j] << endl;
        }
    }
    fout.close();
    fout.open("vorticity.dat");
    for (i = 1; i <= Nx; i ++)
    {
        for (j = 1; j <= Ny; j ++)
        {
            fout << (i - 0.5) * dx << "\t\t" << Ly - (j - 0.5) * dy << "\t\t" << omega[i][j] << endl;
        }
    }
    fout.close();
    fout.open("stream.dat");
    for (i = 1; i <= Nx; i ++)
    {
        for (j = 1; j <= Ny; j ++)
        {
            fout << (i - 0.5) * dx << "\t\t" << Ly - (j - 0.5) * dy << "\t\t" << psi[i][j] << endl;
        }
    }
    fout.close();
    fout.open("u_center.dat");
    for (i = 0; i <= Ny; i ++)
    {
        fout << Ly - i * dy << "\t" << u_center[i] << endl;
    }
    fout.close();
    fout.open("v_center.dat");
    for (i = 0; i <= Nx; i ++)
    {
        fout << i * dx << "\t" << v_center[i] << endl;
    }
    fout.close();
}
int main()
{
    double Re, alpha, uerror_max, verror_max;
    int count = 1, i, j;
    cout << "Enter the Reynolds Number: " ;
    cin >> Re;
    cout << "Enter the relaxation factor for Pressure Poisson Equation: ";
    cin >> alpha;
    initialize(); // Function to initialize velocity and pressure fields
    do
    {
        uerror_max = 0.0;
        verror_max = 0.0;
        cout << "Timestep iteration no: " << count << endl;
        apply_boundary_conditions(); // Apply boundary conditions
        // Solve momentum equations to compute intermediate u and v velocities
        solve_momentum_equations(count, Re);
        solve_pressure_poisson_equation(alpha); // Solve pressure poisson equation using Line Gauss-Siedel meth
        update_velocity(); // Compute corrected velocities
        for (i = 1; i <= Nx; i ++)
        {
            for (j = 1; j <= Ny; j ++)
            {
                if (fabs(u[i][j] - u_old[i][j]) >= uerror_max)
                {
                    uerror_max = fabs(u[i][j] - u_old[i][j]);
                }
                if (fabs(v[i][j] - v_old[i][j]) >= verror_max)
                {
                    verror_max = fabs(v[i][j] - v_old[i][j]);
                }
            }
        }
        if (uerror_max < Tol && verror_max < Tol)
        {
            cout << "Steady state is reached" << endl;
            break;
        }
        for (i = 1; i <= Nx; i ++)
        {
            for (j = 1; j <= Ny; j ++)
            {
                u_old[i][j] = u[i][j];
            }
        }
        for (i = 1; i <= Nx; i ++)
        {
            for (j = 1; j <= Ny; j ++)
            {
                v_old[i][j] = v[i][j];
            }
        }
        count++;
    } while (count <= Max_Timestep);
    compute_velocity_magnitude();
    compute_vorticity();
    compute_stream_function(); // Solve  poisson equation using Line Gauss-Siedel Method
    compute_min_streamline();
    compute_center_velocities();
    outputdata();
    
}