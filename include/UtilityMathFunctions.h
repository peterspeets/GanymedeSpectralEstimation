#ifndef UTILITYMATHFUNCTIONS_H
#define UTILITYMATHFUNCTIONS_H

#include <cmath>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <thread>
#include <sstream>
#include <vector>

#include <kiss_fft.h>
#include <kiss_fftr.h>


template <typename floatingPointType>
class UtilityMathFunctions {

private:


    static uint64_t getTime();

    template <typename T>
    static void saveArrayToFile(const T* array, const int N, const string& filename);
    static void saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename);
    template <typename T>
    static void saveArrayToFile(const complex<T>* cpx, const int N, const string& filename);

    template <typename T>
    struct Spline {
        T a;
        T b;
        T c;
        T d;
        T x0;
        T evaluate(T x) const {
            return a + b*(x - x0)+ c*(x - x0)*(x - x0)+ d*(x - x0)*(x - x0)*(x - x0);
        }
    };



public:
    static complex<floatingPointType>* tvec_gs_i(const floatingPointType* a,const floatingPointType* x,const size_t N,complex<floatingPointType>* y = nullptr);
    static complex<floatingPointType>* tvec_gs_i(const complex<floatingPointType>* a,const floatingPointType* x,const size_t N,complex<floatingPointType>* y = nullptr);
    static complex<floatingPointType>* tvec_gs_i(const complex<floatingPointType>* a,const complex<floatingPointType>* x,const size_t N,complex<floatingPointType>* y = nullptr);
    static floatingPointType** processBScan(floatingPointType** spectra, size_t M,const size_t N, int K,int q_init, int q_i, double vt,int NThreads =1);
    static tuple<complex<floatingPointType>*, floatingPointType> levinson(const complex<floatingPointType>*, size_t N,complex<floatingPointType>* = nullptr );
    static complex<floatingPointType>* polynomialEstimation(const complex<floatingPointType>*, size_t N,complex<floatingPointType>* fa1 = nullptr);
    static pair<floatingPointType*, floatingPointType*> fiaa_oct(const floatingPointType* x, size_t N, int K, int q_i, double vt, floatingPointType* diaaf_floatingPoint = nullptr);
    static void fiaa_oct_partitioned(const floatingPointType* x, size_t N, int K, int numberOfPartitions,int q_i, double vt,floatingPointType* diaaf_floatingPoint = nullptr);
    static void fiaa_oct_loop(floatingPointType** x, int fromIndex, int toIndex,size_t N, int K, int numberOfPartitions,int q_i,
                               double vt,floatingPointType* startingColumn ,floatingPointType** diaaf_floatingPoint);
    static void test();


    class SplineInterpolation {
    public:
        SplineInterpolation(Spline<floatingPointType>** splines , const size_t arraySize );
        floatingPointType evaluate(floatingPointType x);
    //virtual:
        ~SplineInterpolation();
    private:
        const size_t N;
        const Spline<floatingPointType>** splines_ ;
    };


    static SplineInterpolation* splineInterpolation(const floatingPointType*,const floatingPointType*, const size_t);

protected:


};

#include "UtilityMathFunctions.tpp"


#endif // UTILITYMATHFUNCTIONS_H
