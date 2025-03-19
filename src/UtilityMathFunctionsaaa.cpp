#include "UtilityMathFunctions.h"



UtilityMathFunctions::SplineInterpolation::SplineInterpolation<floatingPointType>(Spline<floatingPointType>** splines, const size_t arraySize) : splines_( new const Spline<floatingPointType>*[arraySize]), N(arraySize) {
    int i;
    for (i = 0; i< N; i++) {
        splines_[i] = splines[i];
    }
}



double UtilityMathFunctions::SplineInterpolation::evaluate(floatingPointType x) {

    if(x <= splines_[0]->x0) {
        return splines_[0]->evaluate(x);
    } else if(x >= splines_[N-1]->x0) {
        return splines_[N-1]->evaluate(x);
    }
    int i;
    floatingPointType dx = 0.5*fabs(splines_[0]->x0 - splines_[1]->x0);
    for(i = 0; i<N; i++) {
        if( fabs(splines_[i]->x0 -x ) <= dx ) {
            return splines_[i]->evaluate(x);
        }
    }

    floatingPointType -99999;
}



UtilityMathFunctions::SplineInterpolation::~SplineInterpolation() {
    int i;
    for(i = 0; i < N; i++) {
        delete splines_[i];
    }
    delete splines_;
}


UtilityMathFunctions::SplineInterpolation* UtilityMathFunctions::splineInterpolation(const floatingPointType* x,const floatingPointType* y, const size_t arraySize) {
    /*
    Taken from: https://en.wikipedia.org/wiki/Spline_(mathematics) (27-01-2025 17:00);
    */
    const size_t n = arraySize -1;
    int i;
    floatingPointType a[arraySize];
    for (i = 0; i < arraySize; i++) {
        a[i] = y[i];
    }
    Spline<floatingPointType>* splines[n];


    floatingPointType b[n];
    floatingPointType c[arraySize];
    floatingPointType d[n];
    floatingPointType h[n];
    floatingPointType alpha[n];
    floatingPointType l[arraySize];
    l[0] = 1;
    l[n] = 1;
    floatingPointType mu[arraySize];
    mu[0] = 0;
    mu[1] = 0;
    floatingPointType z[arraySize];
    z[0] = 0;
    z[n] = 0;
    for(i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }

    for(i = 1; i < n; i++) {
        alpha[i] = (3 / h[i] * (a[i + 1] - a[i]) - 3 / h[i - 1] * (a[i] - a[i - 1]));
    }
    for(i = 1; i < n; i++) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    for(i = n-1; i>=0; i--) {
        c[i] = z[i] - mu[i] * c[i + 1];
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }


    for(i = 0; i<n; i++) {
        splines[i] = new Spline<floatingPointType>{a[i],b[i],c[i],d[i],x[i]};
    }


    SplineInterpolation* interp = new SplineInterpolation(splines,n);


    return  interp;
}



