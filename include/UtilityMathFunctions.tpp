
#include "UtilityMathFunctions.h"


template <typename floatingPointType>
template <typename T>
void UtilityMathFunctions<floatingPointType>::saveArrayToFile(const T* array, const int N, const string& filename) {
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << array[i] << endl;
    }
    outFile.close();
}

template <typename floatingPointType>
void UtilityMathFunctions<floatingPointType>::saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename) {
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].r ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile << endl;
    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].i ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile.close();
}




template <typename floatingPointType>
UtilityMathFunctions<floatingPointType>::SplineInterpolation::SplineInterpolation(Spline<floatingPointType>** splines, const size_t arraySize) : splines_( new const Spline<floatingPointType>*[arraySize]), N(arraySize) {
    int i;
    for (i = 0; i< N; i++) {
        splines_[i] = splines[i];
    }
}


template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::tvec_gs_i(const floatingPointType* a,const floatingPointType* x,const size_t N){
    complex<floatingPointType> aComplex[N];
    complex<floatingPointType> xComplex[N];
    for(int i = 0; i < N; i++){
        aComplex[i] = a[i];
        xComplex[i] = x[i];
    }

    return UtilityMathFunctions<floatingPointType>::tvec_gs_i(aComplex,xComplex,N);

}

template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::tvec_gs_i(const complex<floatingPointType>* a,const complex<floatingPointType>* x,const size_t N){
    int i, j;

    kiss_fft_cpx t1t[2*N];
    kiss_fft_cpx t1[2*N];
    //kiss_fft_cpx t[2*N];
    kiss_fft_cpx w[2*N];

    kiss_fft_cpx T1t[2*N];
    kiss_fft_cpx T1[2*N];
    kiss_fft_cpx W[2*N];



    kiss_fft_cpx tmp1[2*N];
    kiss_fft_cpx tmp2[2*N];
    kiss_fft_cpx ftmp[2*N];


    complex<floatingPointType> y[N];
    complex<floatingPointType> z[N];
    //complex<floatingPointType> = t1t[2*N];
    //complex<floatingPointType> = w[2*N];

    kiss_fft_cfg cfg = kiss_fft_alloc(2*N, 0, NULL, NULL);
    kiss_fft_cfg icfg = kiss_fft_alloc(2*N, 1, NULL, NULL);

    for(i = 0; i < N; i++){
        t1[i].r = a[i].real();
        t1[i].i = a[i].imag();
        w[i].r = x[i].real();
        w[i].i = x[i].imag();
    }
    t1[N].r = a[0].real();
    t1[N].i = a[0].imag();
    w[N].r = 0.0;
    w[N].i = 0.0;
    for(i = N+1; i < 2*N; i++){
        t1[i].r = 0.0;
        t1[i].i = 0.0;
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    t1t[0].r = a[0].real();
    t1t[0].i = -a[0].imag();
    for(i = 1; i < N; i++){
        t1t[i].r = 0.0;
        t1t[i].i = 0.0;
    }
    t1t[N].r = a[0].real();
    t1t[N].i = -a[0].imag();
    for(i = N+1; i < 2*N; i++){
        t1t[i].r = a[ 2*N-i ].real();
        t1t[i].i = a[ 2*N-i ].imag();
    }




    kiss_fft( cfg , t1 , T1);
    kiss_fft( cfg , t1t , T1t);
    kiss_fft( cfg , w , W);



    for(i = 0; i < 2*N; i++){
        ftmp[i].r = T1t[i].r*W[i].r - T1t[i].i*W[i].i  ;
        ftmp[i].i = T1t[i].r*W[i].i + T1t[i].i*W[i].r  ;
    }


    kiss_fft( icfg , ftmp , w);



    for(i = N; i < 2*N; i++){
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    kiss_fft( cfg , w , W);
    for(i = 0; i < 2*N; i++){
        ftmp[i].r = T1[i].r*W[i].r - T1[i].i*W[i].i  ;
        ftmp[i].i = T1[i].r*W[i].i + T1[i].i*W[i].r  ;
    }

    kiss_fft( icfg , ftmp, tmp1);




    z[0] = 0.0;
    for(i = 1; i < N;i++){
        z[i] = a[N-i];
    }


    for(i = 0; i < N; i++){
        t1[i].r = z[i].real();
        t1[i].i = z[i].imag();
        w[i].r = x[i].real();
        w[i].i = x[i].imag();
    }
    t1[N].r = z[0].real();
    t1[N].i = z[0].imag();
    w[N].r = 0.0;
    w[N].i = 0.0;
    for(i = N+1; i < 2*N; i++){
        t1[i].r = 0.0;
        t1[i].i = 0.0;
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    t1t[0].r = z[0].real();
    t1t[0].i = -z[0].imag();
    for(i = 1; i < N; i++){
        t1t[i].r = 0.0;
        t1t[i].i = 0.0;
    }
    t1t[N].r = z[0].real();
    t1t[N].i = -z[0].imag();
    for(i = N+1; i < 2*N; i++){
        t1t[i].r = z[ 2*N-i ].real();
        t1t[i].i = z[ 2*N-i ].imag();
    }



    kiss_fft( cfg , t1 , T1);
    kiss_fft( cfg , t1t , T1t);
    kiss_fft( cfg , w , W);



    for(i = 0; i < 2*N; i++){
        ftmp[i].r = T1t[i].r*W[i].r - T1t[i].i*W[i].i  ;
        ftmp[i].i = T1t[i].r*W[i].i + T1t[i].i*W[i].r  ;
    }

    kiss_fft( icfg , ftmp , w);
    for(i = N; i < 2*N; i++){
        w[i].r = 0.0;
        w[i].i = 0.0;
    }



    kiss_fft( cfg , w , W);

    for(i = 0; i < 2*N; i++){
        ftmp[i].r = T1[i].r*W[i].r - T1[i].i*W[i].i  ;
        ftmp[i].i = T1[i].r*W[i].i + T1[i].i*W[i].r  ;
    }



    kiss_fft( icfg , ftmp, tmp2);
    kiss_fft_free(cfg);
    kiss_fft_free(icfg);

    for(i = 0; i < N; i++){
        y[i] =complex<floatingPointType>(tmp1[i].r/(4*N*N) - tmp2[i].r/(4*N*N), tmp1[i].i/(4*N*N) - tmp2[i].i/(4*N*N));
    }



    return y;

}








template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::polynomialEstimation(const complex<floatingPointType>* inputVector, size_t N){
    //Algorithm from: Fast and accurate spectral-estimation axial super-resolution optical coherence tomography, J. de Wit et al. (2021)
    //Matlab code at https://zenodo.org/records/5482794 (J. de Wit)
    floatingPointType M[N];
    complex<floatingPointType>* phi = new complex<floatingPointType>[2*N-1];
    complex<floatingPointType> t1[N];
    complex<floatingPointType> t2[N];
    complex<floatingPointType> s1[N];
    complex<floatingPointType> s2[N];

    kiss_fft_cpx w1[2*N];
    kiss_fft_cpx w2[2*N];
    kiss_fft_cpx t1x[2*N];
    kiss_fft_cpx t2x[2*N];

    kiss_fft_cfg cfg = kiss_fft_alloc(2*N, 0, NULL, NULL);

    for(int i = 0; i < N; i++){
        M[i] = i+1;
        t1[i] = inputVector[i];
    }
    t2[0] = 0.0;



    for(int i = 1; i < N; i++){
        t2[i] = conj(inputVector[N - i]);
    }


    for(int i = 0; i < N; i++){
        s1[i] = conj(t1[N-i-1])*M[i];
        s2[i] = conj(t2[N-i-1])*M[i];
    }


    for(int i = 0; i < N; i++){
        w1[i].r = s1[i].real();
        w1[i].i = s1[i].imag();
        w2[i].r = s2[i].real();
        w2[i].i = s2[i].imag();
        t1x[i].r = t1[i].real();
        t1x[i].i = t1[i].imag();
        t2x[i].r = t2[i].real();
        t2x[i].i = t2[i].imag();
    }
    for(int i = N; i < 2*N; i++){
        w1[i].r = 0.0;
        w2[i].r = 0.0;
        t1x[i].r = 0.0;
        t2x[i].r = 0.0;
        w1[i].i = 0.0;
        w2[i].i = 0.0;
        t1x[i].i = 0.0;
        t2x[i].i = 0.0;
    }
    t1x[N].r = t1[0].real();
    t1x[N].i = t1[0].imag();
    t2x[N].r = t2[0].real();
    t2x[N].i = t2[0].imag();


    kiss_fft_cpx T1X[2*N];
    kiss_fft_cpx T2X[2*N];
    kiss_fft_cpx W1[2*N];
    kiss_fft_cpx W2[2*N];
    kiss_fft_cpx tmp1[2*N];
    kiss_fft_cpx tmp2[2*N];

    kiss_fft_cpx fft1[2*N];
    kiss_fft_cpx fft2[2*N];

    kiss_fft(cfg, t1x, T1X);
    kiss_fft(cfg, t2x, T2X);
    kiss_fft(cfg, w1, W1);
    kiss_fft(cfg, w2, W2);

    kiss_fft_free(cfg);




    complex<floatingPointType> f[N];

    for(int i = 0; i < 2*N; i++){
        fft1[i].r = T1X[i].r * W1[i].r - T1X[i].i * W1[i].i;
        fft1[i].i = T1X[i].r * W1[i].i + T1X[i].i * W1[i].r;
        fft2[i].r = T2X[i].r * W2[i].r - T2X[i].i * W2[i].i;
        fft2[i].i = T2X[i].r * W2[i].i + T2X[i].i * W2[i].r;
    }

    cfg = kiss_fft_alloc(2*N, 1, NULL, NULL);

    kiss_fft(cfg, fft1,tmp1);
    kiss_fft(cfg, fft2,tmp2);


    kiss_fft_free(cfg);

    for(int i = 0; i < N; i++){
        //phi[i] = (tmp1[i].r, -tmp1[i].i) - (tmp2[i].r, -tmp2[i].i);
        phi[i] = complex<floatingPointType>(tmp1[i].r - tmp2[i].r, -tmp1[i].i + tmp2[i].i);
    }
    for(int i = 0; i < N-2; i++){
        phi[N+i] = conj( phi[N-i-2 ] );
    }

    return phi;
}



template <typename floatingPointType>
tuple<complex<floatingPointType>*, floatingPointType, complex<floatingPointType>*> UtilityMathFunctions<floatingPointType>::levinson(const complex<floatingPointType>* inputVector, size_t N){

    //Taken from https://github.com/cokelaer/spectrum/blob/master/src/spectrum/levinson.py

    int i, j,k,kj;
    int khalf;
    complex<floatingPointType> r[N];
    for(i = 0; i<N;i++){
        r[i] = inputVector[i];
    }
    double T0 = r[0].real();

    complex<floatingPointType> T[N-1];
    for(i = 0; i< N-1;i++){
        T[i] = r[i+1];
    }
    int M = N - 1;

    complex<floatingPointType>* A = new complex<floatingPointType>[M];
    complex<floatingPointType>* refer = new complex<floatingPointType>[M];
    //fill(begin(A),end(A),0.0);
    //fill(begin(refer),end(refer),0.0);

    for(i = 0; i < M; i++){
        A[i] = 0;
        refer[i] = 0;
    }

    double P = T0;
    complex<floatingPointType> save;
    complex<floatingPointType> temp;
    for(k = 0; k < M;k++){
        save = T[k];
        if(k == 0){
            temp = -save/P;
        }else{
            for(j = 0; j < k;j++){
                save = save + A[j]*T[k-j-1];
            }
            temp = -save/P;
        }
        P = P*(1.0 - (temp.real()*temp.real() + temp.imag()*temp.imag() ));
        if(P < 0){
            cout << "Singular matrix" << endl;
        }
        A[k] = temp;
        refer[k] = temp;
        if(k == 0){
            continue;
        }
        khalf = (k+1)/2;
        for(j = 0; j < khalf;j++){
            kj = k - j - 1;
            save = A[j];
            A[j] = save + temp*A[kj];
            if(j != kj){
                A[kj] += temp*save;
            }
        }
    }
    return  make_tuple(A,P,refer);
}

template <typename floatingPointType>
floatingPointType UtilityMathFunctions<floatingPointType>::SplineInterpolation::evaluate(floatingPointType x) {

    if(x <= splines_[0]->x0) {
        return splines_[0]->evaluate(x);
    } else if(x >= splines_[N-1]->x0) {
        return splines_[N-1]->evaluate(x);
    }
    int i;
    floatingPointType dx;

    for(i = 0; i<N-1; i++) {
        dx = fabs(splines_[i]->x0 - splines_[i+1]->x0);

        if( fabs(splines_[i]->x0 - x ) <= dx  ) {
            return splines_[i]->evaluate(x);
        }
    }
    cout << "Unexpected x value in interpolation." << endl;
    return -999;
}


template <typename floatingPointType>
UtilityMathFunctions<floatingPointType>::SplineInterpolation::~SplineInterpolation() {
    int i;
    for(i = 0; i < N; i++) {
        delete splines_[i];
    }
    delete splines_;
}

template <typename floatingPointType>
typename UtilityMathFunctions<floatingPointType>::SplineInterpolation* UtilityMathFunctions<floatingPointType>::splineInterpolation(const floatingPointType* x,const floatingPointType* y, const size_t arraySize) {
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



