#include <vector>
#include <iostream>
#include <memory>
#include <functional>
#include <cmath>

using namespace std;

double const GAMMA = 1.4;
double const R = 287;

class StructuredGrid1D
{
public:
    vector<double> data;
    double dx;
    int n;

    StructuredGrid1D(double x0, double dx, int n)
    {
        this->data = vector<double>(n, 0.0);
        for (int i = 0; i < n; i++)
        {
            this->data[i] = x0 + dx * i;
        }
    }
};

class Solution
{
public:
    vector<double> rho;
    vector<double> v;
    vector<double> T;

    Solution(int n)
    {
        this->rho = vector<double>(n, 0.0);
        this->v = vector<double>(n, 0.0);
        this->T = vector<double>(n, 0.0);
    }
};

class Solver
{
    shared_ptr<StructuredGrid1D> grid;
    double CFL;

    vector<double> logA;
    double *logA;

    unique_ptr<Solution> sol;
    double *rho;
    double *v;
    double *T;

    vector<tuple<Solution, double>> solutions;

    Solver(shared_ptr<StructuredGrid1D> _grid, double _CFL, unique_ptr<Solution> initConditions, double rho0, double T0)
    {
        this->grid = _grid;
        this->CFL = CFL;
        this->calculateNozzleArea();

        this->sol = move(initConditions);
        this->rho = this->sol->rho.data();
        this->v = this->sol->v.data();
        this->T = this->sol->T.data();
    }

    void run()
    {
        }

    double calculateMinDt()
    {
        auto dts = vector<double>(this->grid->n, 0.0);

        for (int i = 0; i < this->grid->n; i++)
        {
            dts[i] = this->CFL * this->grid->dx / (sqrt(this->T[i]) + this->v[i]);
        }
        return *min_element(dts.begin(), dts.end());
    }

    void calculateBoundaries()
    {
        int const last = this->grid->n - 1;

        // Inlet
        this->v[0] = 2 * this->v[1] - this->v[2];

        // Outlet
        this->rho[last] = 2 * this->rho[last - 1] - this->rho[last - 2];
        this->v[last] = 2 * this->v[last - 1] - this->v[last - 2];
        this->T[last] = 2 * this->T[last - 1] - this->T[last - 2];
    }

    void calculateNozzleArea()
    {
        auto const &xDomain = this->grid->data;
        for (auto const &x : xDomain)
        {
            this->logA.push_back(log(1 + 2.2 * (x - 1.5) * (x - 1.5)));
        }
    }

    double deriveDensity(int idx, function<double(double const *, int const)> squema)
    {
        double frac1 = this->rho[idx] * squema(this->v, idx);
        double frac2 = this->rho[idx] * this->v[idx] * squema(this->logA.data(), idx);
        double frac3 = this->v[idx] * squema(this->rho, idx);
        return (-1) * (frac1 + frac2 + frac3) / this->grid->dx;
    }

    double deriveVelocity(int idx, function<double(double const *, int const)> squema)
    {
        double frac1 = this->v[idx] * squema(this->v, idx);
        double frac2 = (1 / GAMMA) * (squema(this->T, idx) + this->T[idx] / this->rho[idx] * squema(this->rho, idx));
        return (-1) * (frac1 + frac2) / this->grid->dx;
    }

    double deriveTemperature(int idx, function<double(double const *, int const)> squema)
    {
        double frac1 = this->v[idx] * squema(this->T, idx);
        double frac2 = (GAMMA - 1) * this->T[idx] * (squema(this->v, idx) + this->v[idx] * squema(this->logA.data(), idx));
        return (-1) * (frac1 + frac2) / this->grid->dx;
    }

    double forward(double const *var, int const idx)
    {
        return var[idx + 1] - var[idx];
    }

    double rearward(double const *var, int const idx)
    {
        return var[idx] - var[idx - 1];
    }
};

int main()
{
    auto grid1D = StructuredGrid1D(0.0, 0.1, 31);
    auto vecGrid = grid1D.data;
    for (auto const &point : vecGrid)
    {
        std::cout << point << endl;
    }
    return 0;
}