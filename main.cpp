/*
 * main.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: david thiessen <thiessendg@gmail.com>
 */

/*     Adapted from an example by Glenn Fielder
 *     Simple RK4 integration framework
 *     Copyright (c) 2004, Glenn Fiedler
 *     http://www.gaffer.org/articles
 *
 *     common abbreviations (can be combined):
 *         vert - vertical
 *         horz - horizontal
 *         vel - velocity
 *         acc - acceleration
 *         pos - position
 *         init - initial
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;

constexpr double g = -9.80665; // gravity, the vertical acceleration (m/s/s)
constexpr double re = 6371000.0; // mean Radius of Earth in meters
constexpr double a = 0.0; // horizontal acceleration
constexpr double pi = 3.14159265358979323846; // don't rely on M_PI
constexpr double deg2rad = pi / 180.0;

struct State
{   //state consists of p and p'
	double vertPos; //vertical position, y
	double horzPos; //horizontal position x
	double vertVel; //vertical velocity dy/dt or y'
	double horzVel; //horizontal velocity dx/dt or x'
};

struct Derivative //the derivatives of point p, also can be p'
{   // derivatives of p (p') consists of p' and p''
	double vertVel; // vertical velocity dy/dt or y'
	double horzVel; // horizontal velocity dx/dt or x'
	double vertAcc; // vertical acceleration  dvv/dt or d^2y/dt^2 or y''
	double horzAcc; // horizontal acceleration dvx/dt or d^2x/dt^2 or x''
};

//prototypes
double computeVertAcc(const State &);
double computeHorzAcc(const State &);
const Derivative evaluate(const State &);
const Derivative evaluate(const State &, const double, const Derivative &);
void rk4(State &, const double);

int main(int argc, char * argv[])
{
    double initAlt;
	double initVel;
	double firingAngle;
	double deltaTime;
	double finalTime;

	if (argc != 6)
    {
		//we didn't get command line args, so prompt for them
		cout << "Command line arguments error or not provided." << endl;

		cout << "Enter initial altitude/elevation: " << endl;
        cin >> initAlt;
        while (cin.fail())
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter initial altitude/elevation: " << endl;
            cin >> initAlt;
        }

        cout << "Enter firing angle in degrees (0-90): " << endl;
		cin >> firingAngle;
        while (cin.fail() || firingAngle < 0.0 || firingAngle > 90.0)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter firing angle in degrees (0-90): " << endl;
            cin >> firingAngle;
        }

		cout << "Enter initial velocity (m/s): " << endl;
		cin >> initVel;
        while (cin.fail())
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter initial velocity (m/s): " << endl;
            cin >> initVel;
        }

		cout << "Enter the time step (s) per integration: " << endl;
		cin >> deltaTime;
        while (cin.fail() || deltaTime < 0.)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter the time step (s) per integration: " << endl;
            cin >> deltaTime;
        }

		cout << "Enter final time (s): " << endl;
		cin >> finalTime;
        while (cin.fail() || finalTime < 0.)
        {
            cin.clear();
            cin.ignore(); //skip bad input
            cout << "Error - Invalid input." << endl;
            cout << "Enter final time (s): " << endl;
            cin >> finalTime;
        }
    }
	else
    {
        initAlt = strtod(argv[1], NULL);
		initVel = strtod(argv[2], NULL);
		firingAngle = strtod(argv[3], NULL);
		deltaTime = strtod(argv[4], NULL);
		finalTime = strtod(argv[5], NULL);
	}

	firingAngle *= deg2rad;	//angle radians

	State p; //as in Point P
	p.vertPos = initAlt;
	p.vertVel = initVel * sin(firingAngle);
	p.horzPos = 0.0;
	p.horzVel = initVel * cos(firingAngle);

	double currentTime = 0.0;

	cout << fixed;
	cout << showpoint;

	while (currentTime < finalTime && p.vertPos >= 0.0)
    {
        rk4(p, deltaTime);

		currentTime += deltaTime;

        cout << setprecision(5);
		cout << "t = " << currentTime << endl;
		cout << setprecision(9);
		cout << "\ty = " << p.vertPos << "\ty' = " << p.vertVel << endl;
        cout << "\tx = " << p.horzPos << "\tx' = " << p.horzVel << endl;
	}

	cout << "End of simulation..." << endl;

	return 0;
}

double computeVertAcc(const State & state)
{
	//here, we have a constant gravity
	//F=ma, a = g
	//g varies with alt, g = (re/(re+h))^2
	return g * (re / (re + state.vertPos)) * (re / (re + state.vertPos));
	//gravity/acceleration could be dependent on position and/or velocity?
	//return g - state.y *.000123;
}

double computeHorzAcc(const State & state)
{
	//F=ma, a = F/m = 0
	return a;
}

const Derivative evaluate(const State & p)
{
	Derivative pPrime;
	pPrime.vertVel = p.vertVel;//y' or dy/dt
	pPrime.vertAcc = computeVertAcc(p); //y'' or d^2y/dt^2
	pPrime.horzVel = p.horzVel;//x'
	pPrime.horzAcc = computeHorzAcc(p);//x''
	return pPrime;
}

const Derivative evaluate(const State & initP, const double deltaTime,
		const Derivative & pPrime)
{
    //compute an intermediate state based on initP and its derivative pPrime
	State p;
	//y_n = y_i + v_yi * t;
	p.vertPos = initP.vertPos + pPrime.vertVel * deltaTime;
	//v_yn = v_yi + g * t;
	p.vertVel = initP.vertVel + pPrime.vertAcc * deltaTime;
	//same as above comments
	p.horzPos = initP.horzPos + pPrime.horzVel * deltaTime;
	p.horzVel = initP.horzVel + pPrime.horzAcc * deltaTime;

	Derivative pPrimeNew;
	pPrimeNew.vertVel = p.vertVel;
	pPrimeNew.vertAcc = computeVertAcc(p);
	pPrimeNew.horzVel = p.horzVel;
	pPrimeNew.horzAcc = computeHorzAcc(p);

	return pPrimeNew;
}

void rk4(State &state, const double deltaTime)
{
	//integrate via 4th order runge-kutta
	Derivative k1 = evaluate(state);
	Derivative k2 = evaluate(state, deltaTime/2., k1);
	Derivative k3 = evaluate(state, deltaTime/2., k2);
	Derivative k4 = evaluate(state, deltaTime, k3);

	const double deltaVertVel =
            (k1.vertVel + 2. * (k2.vertVel + k3.vertVel) + k4.vertVel) / 6.;
	const double deltaVertAcc =
            (k1.vertAcc + 2. * (k2.vertAcc + k3.vertAcc) + k4.vertAcc) / 6.;
	const double deltaHorzVel =
            (k1.horzVel + 2. * (k2.horzVel + k3.horzVel) + k4.horzVel) / 6.;
	const double deltaHorzAcc =
            (k1.horzAcc + 2. * (k2.horzAcc + k3.horzAcc) + k4.horzAcc) / 6.;

	//the new state
	state.vertPos = state.vertPos + deltaVertVel * deltaTime;
	state.vertVel = state.vertVel + deltaVertAcc * deltaTime;
	state.horzPos = state.horzPos + deltaHorzVel * deltaTime;
	state.horzVel = state.horzVel + deltaHorzAcc * deltaTime;
}
