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
#include <optional>

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
        T a, b, c, d, x0;
        Spline(T a_, T b_, T c_, T d_, T x0_) : a(a_), b(b_), c(c_), d(d_), x0(x0_) {}//This dummy constructor makes the Spline stuct compatible with vector<>.emplace_back()
        T evaluate(T x) const {
            return a + b*(x - x0) + c*(x - x0)*(x - x0) + d*(x - x0)*(x - x0)*(x - x0);
        }
    };



public:
    static complex<floatingPointType>* gohberg(const floatingPointType* a,const floatingPointType* x,const size_t N,
            complex<floatingPointType>* y = nullptr);
    static complex<floatingPointType>* gohberg(const complex<floatingPointType>* a,const floatingPointType* x,const size_t N,
            complex<floatingPointType>* y = nullptr);
    static complex<floatingPointType>* gohberg(const complex<floatingPointType>* a,const complex<floatingPointType>* x,const size_t N,
            complex<floatingPointType>* y = nullptr);

    static vector<complex<floatingPointType>> gohberg(const vector<floatingPointType>& a,const vector<floatingPointType>& x,
            vector<complex<floatingPointType>>& y = nullopt);
    static vector<complex<floatingPointType>> gohberg(const vector<complex<floatingPointType>>& a,const vector<floatingPointType>& x,
            vector<complex<floatingPointType>>& y = nullopt);
    static vector<complex<floatingPointType>> gohberg(const vector<complex<floatingPointType>>& a,const vector<complex<floatingPointType>>& x,
            vector<complex<floatingPointType>>& y = nullopt);


    static tuple<complex<floatingPointType>*, floatingPointType> levinson(const complex<floatingPointType>*, size_t N,
            complex<floatingPointType>* = nullptr );
    static tuple<vector<complex<floatingPointType>>, floatingPointType> levinson(const vector<complex<floatingPointType>>&,
            vector<complex<floatingPointType>>& = nullopt );
    static complex<floatingPointType>* polynomialEstimation(const complex<floatingPointType>*, size_t N,
            complex<floatingPointType>* phi = nullptr);
    static vector<complex<floatingPointType>> polynomialEstimation(const vector<complex<floatingPointType>>& inputVector,
            vector<complex<floatingPointType>>& phi = nullopt);

    static floatingPointType** processBScan(floatingPointType** spectra, size_t M,const size_t N, int K,int q_init, int q_i, double vt,
                                            int NThreads =1);
    static pair<floatingPointType*, floatingPointType*> fiaa_oct(const floatingPointType* x, size_t N, int K, int q_i, double vt,
            floatingPointType* powerSpectrum = nullptr);
    static void fiaa_oct_partitioned(const floatingPointType* x, size_t N, int K, int numberOfPartitions,int q_i, double vt,
                                     floatingPointType* powerSpectrum = nullptr);
    static void fiaa_oct_loop(floatingPointType** x, int fromIndex, int toIndex,size_t N, int K, int numberOfPartitions,int q_i,
                              double vt,floatingPointType* startingColumn,floatingPointType** powerSpectrum);
    static void test();


    class SplineInterpolation {
    public:
        SplineInterpolation(const vector<floatingPointType>& x, const vector<floatingPointType>& y );
        SplineInterpolation(floatingPointType* x, floatingPointType* y, int arrayLength );
        floatingPointType evaluate(floatingPointType x);
        ~SplineInterpolation();
    private:
        const size_t N;
        vector<Spline<floatingPointType>> splines ;
    };

protected:


};

#include "UtilityMathFunctions.tpp"


#endif // UTILITYMATHFUNCTIONS_H
