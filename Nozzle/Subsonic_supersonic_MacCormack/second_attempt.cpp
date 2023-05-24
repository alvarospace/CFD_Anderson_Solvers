#include <vector>
#include <iostream>
#include <memory>
#include <functional>
#include <cmath>
#include <tuple>

using namespace std;

double const GAMMA = 1.4;
double const R = 287;
double const STEPS = 1400;

double forward_FO_scheme(double const *var, int const idx)
{
    return var[idx + 1] - var[idx];
}

double rearward_FO_scheme(double const *var, int const idx)
{
    return var[idx] - var[idx - 1];
}

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
    double time;
    vector<double> logA;

    unique_ptr<Solution> sol;
    double *rho;
    double *v;
    double *T;
    tuple<int, int> internalRange;

    vector<tuple<unique_ptr<Solution>, double>> solutions;

    Solver(shared_ptr<StructuredGrid1D> _grid, double _CFL, unique_ptr<Solution> initConditions, double rho0, double T0)
    {
        this->time = 0.0;
        this->grid = _grid;
        this->internalRange = {1, (this->grid->n - 1)};
        this->CFL = CFL;
        this->calculateNozzleArea();

        this->updateSolution(move(initConditions));
    }

    void run()
    {
        // TODO: iterate the process and set an exit condition
        double dt = this->calculateMinDt();
        auto tmpSol = make_unique<Solution>(this->grid->n);
        auto nextSol = make_unique<Solution>(this->grid->n);
        auto derivatives = make_unique<Solution>(this->grid->n);
        auto derivatives2 = make_unique<Solution>(this->grid->n);

        // MacCormack predictor-corrector integration (2 steps)
        this->calculateDerivatives(derivatives.get(), &forward_FO_scheme);
        this->integrationStep(tmpSol.get(), this->sol.get(), derivatives.get(), dt);
        this->calculateDerivatives(derivatives2.get(), &rearward_FO_scheme);
        this->integrationStep(nextSol.get(), this->sol.get(), derivatives.get(), dt);
        this->time = this->time + dt;

        this->updateSolution(move(nextSol));
    }

    void calculateAndOverrideDerivativesAverage(Solution *derivatives1, Solution *derivatives2)
    {
        for (int i = get<0>(this->internalRange); i < get<1>(this->internalRange); i++)
        {
            derivatives1->rho[i] = 0.5 * (derivatives1->rho[i] + derivatives2->rho[i]);
            derivatives1->v[i] = 0.5 * (derivatives1->v[i] + derivatives2->v[i]);
            derivatives1->T[i] = 0.5 * (derivatives1->T[i] + derivatives2->T[i]);
        }
    }

    void calculateDerivatives(Solution *derivatives, function<double(double const *, int const)> scheme)
    {
        for (int i = get<0>(this->internalRange); i < get<1>(this->internalRange); i++)
        {
            derivatives->rho[i] = this->deriveDensity(i, scheme);
            derivatives->v[i] = this->deriveVelocity(i, scheme);
            derivatives->T[i] = this->deriveTemperature(i, scheme);
        }
    }

    void integrationStep(Solution *newSol, Solution *currSol, Solution *derivatives, double dt)
    {
        for (int i = get<0>(this->internalRange); i < get<1>(this->internalRange); i++)
        {
            newSol->rho[i] = currSol->rho[i] + derivatives->rho[i] * dt;
            newSol->v[i] = currSol->v[i] + derivatives->v[i] * dt;
            newSol->T[i] = currSol->T[i] + derivatives->T[i] * dt;
        }
    }

    void updateSolution(unique_ptr<Solution> newSol)
    {
        this->solutions.push_back({move(this->sol), this->time});
        this->sol = move(newSol);
        this->rho = this->sol->rho.data();
        this->v = this->sol->v.data();
        this->T = this->sol->T.data();
    }

    double
    calculateMinDt()
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

    double deriveDensity(int idx, function<double(double const *, int const)> scheme)
    {
        double frac1 = this->rho[idx] * scheme(this->v, idx);
        double frac2 = this->rho[idx] * this->v[idx] * scheme(this->logA.data(), idx);
        double frac3 = this->v[idx] * scheme(this->rho, idx);
        return (-1) * (frac1 + frac2 + frac3) / this->grid->dx;
    }

    double deriveVelocity(int idx, function<double(double const *, int const)> scheme)
    {
        double frac1 = this->v[idx] * scheme(this->v, idx);
        double frac2 = (1 / GAMMA) * (scheme(this->T, idx) + this->T[idx] / this->rho[idx] * scheme(this->rho, idx));
        return (-1) * (frac1 + frac2) / this->grid->dx;
    }

    double deriveTemperature(int idx, function<double(double const *, int const)> scheme)
    {
        double frac1 = this->v[idx] * scheme(this->T, idx);
        double frac2 = (GAMMA - 1) * this->T[idx] * (scheme(this->v, idx) + this->v[idx] * scheme(this->logA.data(), idx));
        return (-1) * (frac1 + frac2) / this->grid->dx;
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