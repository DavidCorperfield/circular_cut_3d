#ifndef _CIRCULAR_FITTING_H
#define _CIRCULAR_FITTING_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

using namespace std;

namespace CIRCULAR_FITTING {

	//  define precision by commenting out one of the two lines:

	typedef double reals;       //  defines reals as double (standard for scientific calculations)
								//typedef long double reals;  //  defines reals as long double 

								//   Note: long double is an 80-bit format (more accurate, but more memory demanding and slower)

	typedef long long integers;

	//   next define some frequently used constants:

	const reals One = 1.0, Two = 2.0, Three = 3.0, Four = 4.0, Five = 5.0, Six = 6.0, Ten = 10.0;
	//const reals One=1.0L,Two=2.0L,Three=3.0L,Four=4.0L,Five=5.0L,Six=6.0L,Ten=10.0L;
	const reals Pi = 3.141592653589793238462643383L;
	const reals REAL_MAX = numeric_limits<reals>::max();
	const reals REAL_MIN = numeric_limits<reals>::min();
	const reals REAL_EPSILON = numeric_limits<reals>::epsilon();

	//   next define some frequently used functions:

	template<typename T>
	inline T SQR(T t) { return t*t; }



	class Circle
	{
	public:

		// The fields of a Circle
		reals a, b, r, s, g, Gx, Gy;
		int i, j;

		// constructors
		Circle() { a = 0.; b = 0.; r = 1.; s = 0.; i = 0; j = 0; }
		Circle(reals aa, reals bb, reals rr) { a = aa; b = bb; r = rr; }

		// no destructor we didn't allocate memory by hand.
		void Circle::print(void)
		{
			cout << endl;
			cout << setprecision(10) << "center (" << a << "," << b << ")  radius "
				<< r << "  sigma " << s << "  gradient " << g << "  iter " << i << "  inner " << j << endl;
		}
	};


	// Printing routine
	class Data
	{
	public:
		int n;
		reals *X;		//space is allocated in the constructors
		reals *Y;		//space is allocated in the constructors
		reals meanX, meanY;

		// constructors
		Data() {
			n = 0;
			X = new reals[n];
			Y = new reals[n];
			for (int i = 0; i < n; i++)
			{
				X[i] = 0.;
				Y[i] = 0.;
			}
		}
		Data(int N) {
			n = N;
			X = new reals[n];
			Y = new reals[n];

			for (int i = 0; i < n; i++)
			{
				X[i] = 0.;
				Y[i] = 0.;
			}
		}
		Data(int N, reals dataX[], reals dataY[]) {
			n = N;
			X = new reals[n];
			Y = new reals[n];

			for (int i = 0; i < n; i++)
			{
				X[i] = dataX[i];
				Y[i] = dataY[i];
			}
		}

		// routines
		void means(void) {
			meanX = 0.; meanY = 0.;

			for (int i = 0; i < n; i++)
			{
				meanX += X[i];
				meanY += Y[i];
			}
			meanX /= n;
			meanY /= n;
		}
		void center(void) {
			reals sX = 0., sY = 0.;
			int i;

			for (i = 0; i < n; i++)
			{
				sX += X[i];
				sY += Y[i];
			}
			sX /= n;
			sY /= n;

			for (i = 0; i < n; i++)
			{
				X[i] -= sX;
				Y[i] -= sY;
			}
			meanX = 0.;
			meanY = 0.;
		}
		void scale(void) {
			reals sXX = 0., sYY = 0., scaling;
			int i;

			for (i = 0; i < n; i++)
			{
				sXX += X[i] * X[i];
				sYY += Y[i] * Y[i];
			}
			scaling = sqrt((sXX + sYY) / n / Two);

			for (i = 0; i < n; i++)
			{
				X[i] /= scaling;
				Y[i] /= scaling;
			}
		}
		void print(void) {
			cout << endl << "The data set has " << n << " points with coordinates :" << endl;

			for (int i = 0; i < n - 1; i++) cout << setprecision(7) << "(" << X[i] << "," << Y[i] << "), ";

			cout << "(" << X[n - 1] << "," << Y[n - 1] << ")\n";
		}

		// destructors
		~Data() {
			delete[] X;
			delete[] Y;
		}
	};

	reals Sigma(Data& data, Circle& circle)
	{
		reals sum = 0., dx, dy;

		for (int i = 0; i < data.n; i++)
		{
			dx = data.X[i] - circle.a;
			dy = data.Y[i] - circle.b;
			sum += SQR(sqrt(dx*dx + dy*dy) - circle.r);
		}
		return sqrt(sum / data.n);
	}

	Circle CircleFitByHyper(Data& data)
		/*
		Circle fit to a given set of data points (in 2D)

		This is an algebraic fit based on the journal article

		A. Al-Sharadqah and N. Chernov, "Error analysis for circle fitting algorithms",
		Electronic Journal of Statistics, Vol. 3, pages 886-911, (2009)

		It is an algebraic circle fit with "hyperaccuracy" (with zero essential bias).
		The term "hyperaccuracy" first appeared in papers by Kenichi Kanatani around 2006

		Input:  data     - the class of data (contains the given points):

		data.n   - the number of data points
		data.X[] - the array of X-coordinates
		data.Y[] - the array of Y-coordinates

		Output:
		circle - parameters of the fitting circle:

		circle.a - the X-coordinate of the center of the fitting circle
		circle.b - the Y-coordinate of the center of the fitting circle
		circle.r - the radius of the fitting circle
		circle.s - the root mean square error (the estimate of sigma)
		circle.j - the total number of iterations

		This method combines the Pratt and Taubin fits to eliminate the essential bias.

		It works well whether data points are sampled along an entire circle or
		along a small arc.

		Its statistical accuracy is theoretically higher than that of the Pratt fit
		and Taubin fit, but practically they all return almost identical circles
		(unlike the Kasa fit that may be grossly inaccurate).

		It provides a very good initial guess for a subsequent geometric fit.

		Nikolai Chernov  (September 2012)

		*/
	{
		int i, iter, IterMAX = 99;

		reals Xi, Yi, Zi;
		reals Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
		reals A0, A1, A2, A22;
		reals Dy, xnew, x, ynew, y;
		reals DET, Xcenter, Ycenter;

		Circle circle;

		data.means();   // Compute x- and y- sample means (via a function in the class "data") 

						//     computing moments 

		Mxx = Myy = Mxy = Mxz = Myz = Mzz = 0.;

		for (i = 0; i<data.n; i++)
		{
			Xi = data.X[i] - data.meanX;   //  centered x-coordinates
			Yi = data.Y[i] - data.meanY;   //  centered y-coordinates
			Zi = Xi*Xi + Yi*Yi;

			Mxy += Xi*Yi;
			Mxx += Xi*Xi;
			Myy += Yi*Yi;
			Mxz += Xi*Zi;
			Myz += Yi*Zi;
			Mzz += Zi*Zi;
		}
		Mxx /= data.n;
		Myy /= data.n;
		Mxy /= data.n;
		Mxz /= data.n;
		Myz /= data.n;
		Mzz /= data.n;

		//    computing the coefficients of the characteristic polynomial

		Mz = Mxx + Myy;
		Cov_xy = Mxx*Myy - Mxy*Mxy;
		Var_z = Mzz - Mz*Mz;

		A2 = Four*Cov_xy - Three*Mz*Mz - Mzz;
		A1 = Var_z*Mz + Four*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
		A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
		A22 = A2 + A2;

		//    finding the root of the characteristic polynomial
		//    using Newton's method starting at x=0  
		//     (it is guaranteed to converge to the right root)

		for (x = 0., y = A0, iter = 0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
		{
			Dy = A1 + x*(A22 + 16.*x*x);
			xnew = x - y / Dy;
			if ((xnew == x) || (!isfinite(xnew))) break;
			ynew = A0 + xnew*(A1 + xnew*(A2 + Four*xnew*xnew));
			if (abs(ynew) >= abs(y))  break;
			x = xnew;  y = ynew;
		}

		//    computing paramters of the fitting circle

		DET = x*x - x*Mz + Cov_xy;
		Xcenter = (Mxz*(Myy - x) - Myz*Mxy) / DET / Two;
		Ycenter = (Myz*(Mxx - x) - Mxz*Mxy) / DET / Two;

		//       assembling the output

		circle.a = Xcenter + data.meanX;
		circle.b = Ycenter + data.meanY;
		circle.r = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz - x - x);
		circle.s = Sigma(data, circle);
		circle.i = 0;
		circle.j = iter;  //  return the number of iterations, too

		return circle;
	};

}
#endif // _CIRCULAR_FITTING_H