#ifndef OPENMM_SPLINEFITTER_H_
#define OPENMM_SPLINEFITTER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "windowsExport.h"
#include <vector>

namespace OpenMM {

/**
 * SplineFitter provides routines for performing cubic spline interpolation.
 */

class OPENMM_EXPORT SplineFitter {
public:
    /**
     * Fit a cubic spline to a set of data points.  The resulting spline interpolates all the
     * data points and has a continuous second derivative everywhere. The second derivatives are
     * identical at the end points if periodic=true or 0 at the end points if periodic=false.
     *
     * @param x        the values of the independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: x[i] > x[i-1].
     * @param y        the values of the dependent variable at the data points to interpolate
     * @param periodic whether the interpolated function is periodic
     * @param deriv    on exit, this contains the second derivative of the spline at each of the data points
     */
    static void createSpline(const std::vector<double>& x, const std::vector<double>& y, bool periodic, std::vector<double>& deriv);
    /**
     * Fit a natural cubic spline to a set of data points.  The resulting spline interpolates all the
     * data points, has a continuous second derivative everywhere, and has a second derivative of 0 at
     * its end points.
     *
     * @param x      the values of the independent variable at the data points to interpolate.  They must
     *               be strictly increasing: x[i] > x[i-1].
     * @param y      the values of the dependent variable at the data points to interpolate
     * @param deriv  on exit, this contains the second derivative of the spline at each of the data points
     */
    static void createNaturalSpline(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& deriv);
    /**
     * Fit a periodic cubic spline to a set of data points.  The resulting spline interpolates all the
     * data points, has a continuous second derivative everywhere, and has identical second derivatives
     * at the end points.
     *
     * @param x      the values of the independent variable at the data points to interpolate.  They must
     *               be strictly increasing: x[i] > x[i-1].
     * @param y      the values of the dependent variable at the data points to interpolate.  The first and
     *               last entries must be identical.
     * @param deriv  on exit, this contains the second derivative of the spline at each of the data points
     */
    static void createPeriodicSpline(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& deriv);
    /**
     * Evaluate a 1D spline generated by one of the other methods in this class.
     *
     * @param x     the values of the independent variable at the data points to interpolate
     * @param y     the values of the dependent variable at the data points to interpolate
     * @param deriv the vector of second derivatives that was calculated by one of the other methods
     * @param t     the value of the independent variable at which to evaluate the spline
     * @return the value of the spline at the specified point
     */
    static double evaluateSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& deriv, double t);
    /**
     * Evaluate the derivative of a 1D spline generated by one of the other methods in this class.
     *
     * @param x     the values of the independent variable at the data points to interpolate
     * @param y     the values of the dependent variable at the data points to interpolate
     * @param deriv the vector of second derivatives that was calculated by one of the other methods
     * @param t     the value of the independent variable at which to evaluate the spline
     * @return the value of the spline's derivative  at the specified point
     */
    static double evaluateSplineDerivative(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& deriv, double t);
    /**
     * Fit a cubic spline surface f(x,y) to a 2D set of data points.  The resulting spline interpolates all the
     * data points and has a continuous second derivative everywhere. The second derivatives are identical at
     * the boundary if periodic=true or 0 at the boundary if periodic=false.
     *
     * @param x        the values of the first independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: x[i] > x[i-1].
     * @param y        the values of the second independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: y[i] > y[i-1].
     * @param values   the values of the dependent variable at the data points to interpolate.  They must be ordered
     *                 so that values[i+xsize*j] = f(x[i],y[j]), where xsize is the length of x.
     * @param periodic whether the interpolated function is periodic
     * @param c        on exit, this contains the spline coefficients at each of the data points
     */
    static void create2DSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& values, bool periodic, std::vector<std::vector<double> >& c);
    /**
     * Fit a natural cubic spline surface f(x,y) to a 2D set of data points.  The resulting spline interpolates all the
     * data points, has a continuous second derivative everywhere, and has a second derivative of 0 at the boundary.
     *
     * @param x      the values of the first independent variable at the data points to interpolate.  They must
     *               be strictly increasing: x[i] > x[i-1].
     * @param y      the values of the second independent variable at the data points to interpolate.  They must
     *               be strictly increasing: y[i] > y[i-1].
     * @param values the values of the dependent variable at the data points to interpolate.  They must be ordered
     *               so that values[i+xsize*j] = f(x[i],y[j]), where xsize is the length of x.
     * @param c      on exit, this contains the spline coefficients at each of the data points
     */
    static void create2DNaturalSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& values, std::vector<std::vector<double> >& c);
    /**
     * Evaluate a 2D spline generated by one of the other methods in this class.
     *
     * @param x      the values of the first independent variable at the data points to interpolate
     * @param y      the values of the second independent variable at the data points to interpolate
     * @param values the values of the dependent variable at the data points to interpolate
     * @param c      the vector of spline coefficients that was calculated by one of the other methods
     * @param u      the value of the first independent variable at which to evaluate the spline
     * @param v      the value of the second independent variable at which to evaluate the spline
     * @return the value of the spline at the specified point
     */
    static double evaluate2DSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& values, const std::vector<std::vector<double> >& c, double u, double v);
    /**
     * Evaluate the derivatives of a 2D spline generated by one of the other methods in this class.
     *
     * @param x      the values of the first independent variable at the data points to interpolate
     * @param y      the values of the second independent variable at the data points to interpolate
     * @param values the values of the dependent variable at the data points to interpolate
     * @param c      the vector of spline coefficients that was calculated by one of the other methods
     * @param u      the value of the first independent variable at which to evaluate the spline
     * @param v      the value of the second independent variable at which to evaluate the spline
     * @param dx     on exit, the x derivative of the spline at the specified point
     * @param dy     on exit, the y derivative of the spline at the specified point
     */
    static void evaluate2DSplineDerivatives(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& values, const std::vector<std::vector<double> >& c, double u, double v, double& dx, double& dy);
    /**
     * Fit a cubic spline surface f(x,y,z) to a 3D set of data points.  The resulting spline interpolates all the
     * data points and has a continuous second derivative everywhere. The second derivatives are identical at
     * the boundary if periodic=true or 0 at the boundary if periodic=false.
     *
     * @param x        the values of the first independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: x[i] > x[i-1].
     * @param y        the values of the second independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: y[i] > y[i-1].
     * @param z        the values of the third independent variable at the data points to interpolate.  They must
     *                 be strictly increasing: z[i] > z[i-1].
     * @param values   the values of the dependent variable at the data points to interpolate.  They must be ordered
     *                 so that values[i+xsize*j+xsize*ysize*k] = f(x[i],y[j],z[k]), where xsize is the length of x
     *                 and ysize is the length of y.
     * @param periodic whether the interpolated function is periodic
     * @param c        on exit, this contains the spline coefficients at each of the data points
     */
    static void create3DSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& values, bool periodic, std::vector<std::vector<double> >& c);
    /**
     * Fit a natural cubic spline surface f(x,y,z) to a 3D set of data points.  The resulting spline interpolates all the
     * data points, has a continuous second derivative everywhere, and has a second derivative of 0 at the boundary.
     *
     * @param x      the values of the first independent variable at the data points to interpolate.  They must
     *               be strictly increasing: x[i] > x[i-1].
     * @param y      the values of the second independent variable at the data points to interpolate.  They must
     *               be strictly increasing: y[i] > y[i-1].
     * @param z      the values of the third independent variable at the data points to interpolate.  They must
     *               be strictly increasing: z[i] > z[i-1].
     * @param values the values of the dependent variable at the data points to interpolate.  They must be ordered
     *               so that values[i+xsize*j+xsize*ysize*k] = f(x[i],y[j],z[k]), where xsize is the length of x
     *               and ysize is the length of y.
     * @param c      on exit, this contains the spline coefficients at each of the data points
     */
    static void create3DNaturalSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& values, std::vector<std::vector<double> >& c);
    /**
     * Evaluate a 3D spline generated by one of the other methods in this class.
     *
     * @param x      the values of the first independent variable at the data points to interpolate
     * @param y      the values of the second independent variable at the data points to interpolate
     * @param z      the values of the third independent variable at the data points to interpolate
     * @param values the values of the dependent variable at the data points to interpolate
     * @param c      the vector of spline coefficients that was calculated by one of the other methods
     * @param u      the value of the first independent variable at which to evaluate the spline
     * @param v      the value of the second independent variable at which to evaluate the spline
     * @param w      the value of the third independent variable at which to evaluate the spline
     * @return the value of the spline at the specified point
     */
    static double evaluate3DSpline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& values, const std::vector<std::vector<double> >& c, double u, double v, double w);
    /**
     * Evaluate the derivatives of a 3D spline generated by one of the other methods in this class.
     *
     * @param x      the values of the first independent variable at the data points to interpolate
     * @param y      the values of the second independent variable at the data points to interpolate
     * @param z      the values of the third independent variable at the data points to interpolate
     * @param values the values of the dependent variable at the data points to interpolate
     * @param c      the vector of spline coefficients that was calculated by one of the other methods
     * @param u      the value of the first independent variable at which to evaluate the spline
     * @param v      the value of the second independent variable at which to evaluate the spline
     * @param w      the value of the third independent variable at which to evaluate the spline
     * @param dx     on exit, the x derivative of the spline at the specified point
     * @param dy     on exit, the y derivative of the spline at the specified point
     * @param dz     on exit, the z derivative of the spline at the specified point
     */
    static void evaluate3DSplineDerivatives(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& values, const std::vector<std::vector<double> >& c, double u, double v, double w, double& dx, double& dy, double &dz);
private:
    static void solveTridiagonalMatrix(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs, std::vector<double>& sol);
};

} // namespace OpenMM

#endif /*OPENMM_SPLINEFITTER_H_*/
