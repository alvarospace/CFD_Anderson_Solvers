#include <vector>
#include <cmath>
#include <functional>
#include <tuple>
#include <algorithm>

using namespace std;

double const GAMMA = 1.4;
double const R = 287;

double calculareNozzleArea(double const x)
{
    return 1 + 2.2 * (x - 1.5) * (x - 1.5);
}

class Point
{
public:
    double rho;
    double v;
    double T;
    double x;
    double A;
    double a;

    Point(double _density, double _velocity, double _temperature, double _x)
    {
        this->rho = _density;
        this->v = _velocity;
        this->T = _temperature;
        this->x = _x;
        this->A = this->calculateNozzleArea(_x);
        this->a = this->calculateSoundVelocity(_temperature);
    }

private:
    double calculateNozzleArea(double const x)
    {
        return 1 + 2.2 * (x - 1.5) * (x - 1.5);
    }

    double calculateSoundVelocity(double T)
    {
        return sqrt(GAMMA * R * T);
    }
};

class PointDerivatives
{
public:
    double dRho;
    double dV;
    double dT;

    PointDerivatives(double _dRho, double _dV, double _dT)
    {
        this->dRho = _dRho;
        this->dV = _dV;
        this->dT = _dT;
    }
};

class Solution
{
public:
    vector<Point> points;

    Solution(int const numPoints)
    {
        this->points = vector<Point>(numPoints, Point(0, 0, 0, 0));
    }
};

class Solver
{
    double dx;
    int gridPoints;
    vector<double> grid;
    double courant;
    vector<tuple<Solution, double>> solutions;

    double const xLim = 3.0;

    Solver(double deltaX, double _courant)
    {
        this->dx = deltaX;
        this->gridPoints = floor(xLim / dx);
        this->courant = _courant;
    }

    void setInitConditions(Solution initCond)
    {
        this->solutions.push_back(make_tuple(initCond, 0.0));
    }

    void run()
    {
        // TODO: boundaries, convergence criteria, MacCormack 2 steps
        Solution initCond = get<0>(this->solutions[0]);
    }

    void advanceSolution(Solution &sol)
    {
        double dt = this->calculateMinDt(sol);
    }

    double calculateMinDt(Solution &sol)
    {
        vector<double> dts;

        for (auto const &point : sol.points)
        {
            dts.push_back(this->courant * this->dx / (point.a + point.v));
        }
        return *min_element(dts.begin(), dts.end());
    }

    Point integratePoint(double dt, Point &target, PointDerivatives &derivatives)
    {
        double rho = target.rho + derivatives.dRho * dt;
        double v = target.v + derivatives.dV * dt;
        double T = target.T + derivatives.dT * dt;

        return Point(rho, v, T, target.x);
    }

    double deriveDensity(Point &local, Point &next, function<double(double const &, double const &)> squema)
    {
        double frac1 = local.rho * squema(local.v, next.v);
        double frac2 = local.rho * local.v * squema(log(local.A), log(next.A));
        double frac3 = local.v * squema(local.rho, next.rho);
        return (-1) * (frac1 + frac2 + frac3) / this->dx;
    }

    double deriveVelocity(Point &local, Point &next, function<double(double const &, double const &)> squema)
    {
        double frac1 = local.v * squema(local.v, next.v);
        double frac2 = (1 / GAMMA) * (squema(local.T, next.T) + local.T / local.rho * squema(local.rho, next.rho));
        return (-1) * (frac1 + frac2) / this->dx;
    }

    double deriveTemperature(Point &local, Point &next, function<double(double const &, double const &)> squema)
    {
        double frac1 = local.v * squema(local.T, next.T);
        double frac2 = (GAMMA - 1) * local.T * (squema(local.v, next.v) + local.v * squema(log(local.A), log(next.A)));
        return (-1) * (frac1 + frac2) / this->dx;
    }

    double forward(double const &local, double const &next)
    {
        return next - local;
    }

    double rearward(double const &local, double const &prev)
    {
        return local - prev;
    }
};
