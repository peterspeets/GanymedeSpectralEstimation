
#include "UtilityMathFunctions.h"
template <typename floatingPointType>
uint64_t UtilityMathFunctions<floatingPointType>::getTime() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

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
template <typename T>
void UtilityMathFunctions<floatingPointType>::saveArrayToFile(const complex<T>* cpx, const int N, const string& filename) {
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].real() ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile << endl;
    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].imag() ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile.close();

}



template <typename floatingPointType>
floatingPointType** UtilityMathFunctions<floatingPointType>::processBScan(floatingPointType** spectra, size_t M,const size_t N, int K, int q_i, double vt) {
    int i, j;
    floatingPointType** processedImage = new floatingPointType*[M];
    //cout << "Performing FIAA: " << endl;

    uint64_t time0;
    uint64_t time1;

    pair<floatingPointType*, floatingPointType*> fiaa_output;
    for(i = 0; i < M; i++) {

        time0 = getTime();
        fiaa_output = fiaa_oct(spectra[i],N,K,q_i,vt);
        time1  = getTime();
        //cout << "fiaa time: " << (time1 -time0)/q_i << endl;

        cout << i << endl;
        processedImage[i] = fiaa_output.first;

    }

    return processedImage;

}


template <typename floatingPointType>
pair<floatingPointType*, floatingPointType*> UtilityMathFunctions<floatingPointType>::fiaa_oct(const floatingPointType* x, size_t N, int K, int q_i, double vt) {
    int i, j;

    ostringstream  filename;

    floatingPointType* Eta = new floatingPointType[q_i+1];
    floatingPointType eta = 0.0;
    floatingPointType af;
    for(i = 0; i < N; i++) {
        eta += abs(x[i]*x[i]);
    }
    eta /= N;


    Eta[0] = eta;
    kiss_fft_cpx diaaf[K];
    floatingPointType* diaaf_floatingPoint = new floatingPointType[K];
    kiss_fft_cpx temp[K];
    kiss_fft_cpx temp2[K];
    kiss_fft_cpx q[K];
    kiss_fft_cpx Fa1[K];

    complex<floatingPointType> c[N];
    kiss_fft_cpx diaa_num[K];
    complex<floatingPointType> diaa_den[K];
    floatingPointType diag_a[N];
    complex<floatingPointType>* A;
    complex<floatingPointType>* y;
    complex<floatingPointType>* fa1;

    uint64_t time0;
    uint64_t time1;
    uint64_t time3;
    uint64_t time4;


    kiss_fft_cfg cfg = kiss_fft_alloc(K, 0, NULL, NULL);
    kiss_fft_cfg icfg = kiss_fft_alloc(K, 1, NULL, NULL);


    for(i = 0; i<N; i++) {
        temp[i].r = x[i];
        temp[i].i = 0.0;
    }
    for(i = N; i<K; i++) {
        temp[i].r = 0.0;
        temp[i].i = 0.0;
    }

    kiss_fft( cfg, temp, diaaf);
    for(i = 0; i<K; i++) {
        diaaf[i].r = (diaaf[i].r *diaaf[i].r  + diaaf[i].i *diaaf[i].i)/(N*N);
        diaaf[i].i=0;
    }




    for(int k = 0; k < q_i; k++) {
        kiss_fft(icfg,diaaf,q);

        //filename.str("");
        //filename << "D:\\data\\q" << k << ".txt";
        //saveArrayToFile(q, K, filename.str());

        for(i = 0; i<N; i++) {
            c[i] = complex<floatingPointType>(q[i].r, 0.0*q[i].i);
        }
        c[0] += vt*eta;

        //cout << "eta: "<< eta << endl;

        time0 = getTime();

        tuple<complex<floatingPointType>*, floatingPointType> levinsonOut = levinsonUnsafe(c,N);
        time1  = getTime();
        //cout << "levinson time: " << (time1 -time0) << endl;


        //delete[] get<2>(levinsonOut);
        A = get<0>(levinsonOut);
        af = sqrt(get<1>(levinsonOut));

        for(i = 0; i < N; i++) {
            A[i] /= af;
        }
        time0  = getTime();
        y = tvec_gs_i(A,x,N);
        time1  = getTime();
        //cout << "tvec time: " << (time1 -time0) << endl;

        for(i = 0; i<N; i++) {
            temp[i].r = y[i].real();
            temp[i].i = y[i].imag();
        }
        for(i =  N; i<K; i++) {
            temp[i].r = 0;
            temp[i].i = 0;
        }
        kiss_fft(cfg,temp,diaa_num);




        time0 = getTime();
        fa1 = polynomialEstimationLowMem(A,N);//fa1 has size 2*N-1




        time1  = getTime();
        //cout << "poly time: " << (time1 -time0) << endl;


        for(i = 0; i < 2*N-1; i++) {
            temp[i].r = fa1[i].real();
            temp[i].i = fa1[i].imag();
        }

        kiss_fft(cfg,temp,Fa1);


        diaa_den[0] = complex<floatingPointType>(Fa1[0].r, Fa1[0].i);

        for(i = 1; i<K; i++) {
            diaa_den[i] = complex<floatingPointType>(Fa1[ K - i ].r, Fa1[K - i  ].i);
        }

        for(i = 0; i<K; i++) {
            diaaf[i].r  = abs( (diaa_num[i].r*diaa_num[i].r + diaa_num[i].i*diaa_num[i].i  )  /(diaa_den[i]*conj(diaa_den[i])) );
            diaaf[i].i = 0.0;
        }


        for(i = 1; i < N; i++) {
            temp[i].r = A[N-i].real();
            temp[i].i = -A[N-i].imag();
        }
        temp[0].r = 0.0;
        temp[0].i = 0.0;

        //saveArrayToFile(diaa_num, K, "D:\\data\\num.txt");
        //saveArrayToFile(diaa_den, K, "D:\\data\\den.txt");


        diag_a[0] = A[0].real()*A[0].real() + A[0].imag()*A[0].imag()  -(temp[0].r * temp[0].r +  temp[0].i*temp[0].i);

        for(i = 1; i < N; i++) {
            diag_a[i] = diag_a[i-1] +  A[i].real()*A[i].real() + A[i].imag()*A[i].imag()  -(temp[i].r * temp[i].r +  temp[i].i*temp[i].i);
        }


        eta = 0.0;
        for(i = 0; i < N; i++) {
            eta +=  abs(y[i]*y[i] / (diag_a[i]*diag_a[i])) ;
        }
        eta /= N;
        Eta[k] = eta;
    }


    kiss_fft_free(cfg);
    kiss_fft_free(icfg);
    delete[] A;
    delete[] y;
    delete[] fa1;
    pair<floatingPointType*, floatingPointType*> result;

    for(i = 0; i < K; i++) {
        diaaf_floatingPoint[i] = diaaf[i].r;
    }

    result.first = diaaf_floatingPoint;
    result.second = Eta;
    return result;
}


template <typename floatingPointType>
tuple<complex<floatingPointType>*, floatingPointType> UtilityMathFunctions<floatingPointType>::levinsonUnsafe(const complex<floatingPointType>* inputVector, size_t N){
    int i, j,k,kj;
    int khalf;
    complex<floatingPointType> r[N];

    floatingPointType T0 = inputVector[0].real();

    int M = N - 1;

    complex<floatingPointType>* A = new complex<floatingPointType>[N];

    A[0] = 1.0;
    for(i = 0; i < M; i++) {
        A[i+1] = 0;
    }

    floatingPointType P = T0;
    complex<floatingPointType> save;
    complex<floatingPointType> temp;

    bool warnedSingularMatrix = false;

    for(k = 0; k < M; k++) {
        save = inputVector[k+1];
        if(k == 0) {
            temp = -save/P;
        } else {
            for(j = 0; j < k; j++) {
                save = save + A[j+1]*inputVector[k-j];
            }
            temp = -save/P;
        }
        P = P*(1.0 - (temp.real()*temp.real() + temp.imag()*temp.imag() ));
        if(P < 0 && !warnedSingularMatrix) {
            warnedSingularMatrix = true;
            cout << "Singular matrix " << endl;

            saveArrayToFile(inputVector, N, "D:/data/OffendingVector.txt");


        }
        A[k+1] = temp;
        if(k == 0) {
            continue;
        }
        khalf = (k+1)/2;
        for(j = 0; j < khalf; j++) {
            kj = k - j - 1;
            save = A[j+1];
            A[j+1] = save + temp*A[kj+1];
            if(j != kj) {
                A[kj+1] += temp*save;
            }
        }
    }
    return  make_tuple(A,P);
}



template <typename floatingPointType>
UtilityMathFunctions<floatingPointType>::SplineInterpolation::SplineInterpolation(Spline<floatingPointType>** splines, const size_t arraySize) : splines_( new const Spline<floatingPointType>*[arraySize]), N(arraySize) {
    int i;
    for (i = 0; i< N; i++) {
        splines_[i] = splines[i];
    }
}


template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::tvec_gs_i(const floatingPointType* a,const floatingPointType* x,const size_t N) {
    complex<floatingPointType> aComplex[N];
    complex<floatingPointType> xComplex[N];
    for(int i = 0; i < N; i++) {
        aComplex[i] = a[i];
        xComplex[i] = x[i];
    }

    return UtilityMathFunctions<floatingPointType>::tvec_gs_i(aComplex,xComplex,N);
}

template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::tvec_gs_i(const complex<floatingPointType>* a,const floatingPointType* x,const size_t N) {
    complex<floatingPointType> xComplex[N];
    for(int i = 0; i < N; i++) {
        xComplex[i] = x[i];
    }
    return UtilityMathFunctions<floatingPointType>::tvec_gs_i(a,xComplex,N);
}


template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::tvec_gs_i(const complex<floatingPointType>* a,const complex<floatingPointType>* x,const size_t N) {
    int i, j;


    kiss_fft_cpx t1t[2*N];
    kiss_fft_cpx t1[2*N];
    kiss_fft_cpx w[2*N];

    kiss_fft_cpx T1t[2*N];
    kiss_fft_cpx T1[2*N];
    kiss_fft_cpx W[2*N];



    kiss_fft_cpx tmp1[2*N];
    kiss_fft_cpx tmp2[2*N];
    kiss_fft_cpx ftmp[2*N];


    complex<floatingPointType>* y = new complex<floatingPointType>[N];
    complex<floatingPointType> z[N];

    kiss_fft_cfg cfg = kiss_fft_alloc(2*N, 0, NULL, NULL);
    kiss_fft_cfg icfg = kiss_fft_alloc(2*N, 1, NULL, NULL);

    for(i = 0; i < N; i++) {
        t1[i].r = a[i].real();
        t1[i].i = a[i].imag();
        w[i].r = x[i].real();
        w[i].i = x[i].imag();
    }
    t1[N].r = a[0].real();
    t1[N].i = a[0].imag();
    w[N].r = 0.0;
    w[N].i = 0.0;
    for(i = N+1; i < 2*N; i++) {
        t1[i].r = 0.0;
        t1[i].i = 0.0;
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    t1t[0].r = a[0].real();
    t1t[0].i = -a[0].imag();
    for(i = 1; i < N; i++) {
        t1t[i].r = 0.0;
        t1t[i].i = 0.0;
    }
    t1t[N].r = a[0].real();
    t1t[N].i = -a[0].imag();
    for(i = N+1; i < 2*N; i++) {
        t1t[i].r = a[ 2*N-i ].real();
        t1t[i].i = a[ 2*N-i ].imag();
    }


    kiss_fft( cfg, t1, T1);
    kiss_fft( cfg, t1t, T1t);
    kiss_fft( cfg, w, W);

    //saveArrayToFile(t1, 2*N,"D:\\data\\testt1.txt");
    //saveArrayToFile(t1t, 2*N,"D:\\data\\testt1t.txt");
    //saveArrayToFile(w, 2*N,"D:\\data\\testw.txt");

    for(i = 0; i < 2*N; i++) {
        ftmp[i].r = T1t[i].r*W[i].r - T1t[i].i*W[i].i  ;
        ftmp[i].i = T1t[i].r*W[i].i + T1t[i].i*W[i].r  ;
    }
    kiss_fft( icfg, ftmp, w);



    for(i = N; i < 2*N; i++) {
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    kiss_fft( cfg, w, W);
    for(i = 0; i < 2*N; i++) {
        ftmp[i].r = T1[i].r*W[i].r - T1[i].i*W[i].i  ;
        ftmp[i].i = T1[i].r*W[i].i + T1[i].i*W[i].r  ;
    }

    kiss_fft( icfg, ftmp, tmp1);




    z[0] = 0.0;
    for(i = 1; i < N; i++) {
        z[i] = a[N-i];
    }


    for(i = 0; i < N; i++) {
        t1[i].r = z[i].real();
        t1[i].i = z[i].imag();
        w[i].r = x[i].real();
        w[i].i = x[i].imag();
    }
    t1[N].r = z[0].real();
    t1[N].i = z[0].imag();
    w[N].r = 0.0;
    w[N].i = 0.0;
    for(i = N+1; i < 2*N; i++) {
        t1[i].r = 0.0;
        t1[i].i = 0.0;
        w[i].r = 0.0;
        w[i].i = 0.0;
    }

    t1t[0].r = z[0].real();
    t1t[0].i = -z[0].imag();
    for(i = 1; i < N; i++) {
        t1t[i].r = 0.0;
        t1t[i].i = 0.0;
    }
    t1t[N].r = z[0].real();
    t1t[N].i = -z[0].imag();
    for(i = N+1; i < 2*N; i++) {
        t1t[i].r = z[ 2*N-i ].real();
        t1t[i].i = z[ 2*N-i ].imag();
    }



    kiss_fft( cfg, t1, T1);
    kiss_fft( cfg, t1t, T1t);
    kiss_fft( cfg, w, W);



    for(i = 0; i < 2*N; i++) {
        ftmp[i].r = T1t[i].r*W[i].r - T1t[i].i*W[i].i  ;
        ftmp[i].i = T1t[i].r*W[i].i + T1t[i].i*W[i].r  ;
    }

    kiss_fft( icfg, ftmp, w);
    for(i = N; i < 2*N; i++) {
        w[i].r = 0.0;
        w[i].i = 0.0;
    }



    kiss_fft( cfg, w, W);

    for(i = 0; i < 2*N; i++) {
        ftmp[i].r = T1[i].r*W[i].r - T1[i].i*W[i].i  ;
        ftmp[i].i = T1[i].r*W[i].i + T1[i].i*W[i].r  ;
    }


    kiss_fft( icfg, ftmp, tmp2);
    kiss_fft_free(cfg);
    kiss_fft_free(icfg);

    for(i = 0; i < N; i++) {
        y[i] =complex<floatingPointType>(tmp1[i].r/(4*N*N) - tmp2[i].r/(4*N*N), tmp1[i].i/(4*N*N) - tmp2[i].i/(4*N*N));
    }



    return y;

}








template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::polynomialEstimation(const complex<floatingPointType>* inputVector, size_t N) {
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

    for(int i = 0; i < N; i++) {
        M[i] = i+1;
        t1[i] = inputVector[i];
    }
    t2[0] = 0.0;



    for(int i = 1; i < N; i++) {
        t2[i] = conj(inputVector[N - i]);
    }


    for(int i = 0; i < N; i++) {
        s1[i] = conj(t1[N-i-1])*M[i];
        s2[i] = conj(t2[N-i-1])*M[i];
    }


    for(int i = 0; i < N; i++) {
        w1[i].r = s1[i].real();
        w1[i].i = s1[i].imag();
        w2[i].r = s2[i].real();
        w2[i].i = s2[i].imag();
        t1x[i].r = t1[i].real();
        t1x[i].i = t1[i].imag();
        t2x[i].r = t2[i].real();
        t2x[i].i = t2[i].imag();
    }


    for(int i = N; i < 2*N; i++) {
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

    for(int i = 0; i < 2*N; i++) {
        fft1[i].r = T1X[i].r * W1[i].r - T1X[i].i * W1[i].i;
        fft1[i].i = T1X[i].r * W1[i].i + T1X[i].i * W1[i].r;
        fft2[i].r = T2X[i].r * W2[i].r - T2X[i].i * W2[i].i;
        fft2[i].i = T2X[i].r * W2[i].i + T2X[i].i * W2[i].r;
    }

    cfg = kiss_fft_alloc(2*N, 1, NULL, NULL);

    kiss_fft(cfg, fft1,tmp1);
    kiss_fft(cfg, fft2,tmp2);


    kiss_fft_free(cfg);

    for(int i = 0; i < N; i++) {
        phi[i] = complex<floatingPointType>( (tmp1[i].r - tmp2[i].r)/(2*N), ( -tmp1[i].i + tmp2[i].i)/(2*N)    );
    }
    for(int i = 0; i < N-2; i++) {
        phi[N+i] = conj( phi[N-i-2 ] );
    }

    return phi;
}






template <typename floatingPointType>
complex<floatingPointType>* UtilityMathFunctions<floatingPointType>::polynomialEstimationLowMem(const complex<floatingPointType>* inputVector, size_t N) {
    //Algorithm from: Fast and accurate spectral-estimation axial super-resolution optical coherence tomography, J. de Wit et al. (2021)
    //Matlab code at https://zenodo.org/records/5482794 (J. de Wit)

    complex<floatingPointType> t[N];
    complex<floatingPointType> s[N];


    kiss_fft_cpx w[2*N];
    kiss_fft_cpx tx[2*N];

    kiss_fft_cfg cfg = kiss_fft_alloc(2*N, 0, NULL, NULL);

    for(int i = 0; i < N; i++) {
        t[i] = inputVector[i];
    }


    for(int i = 0; i < N; i++) {
        s[i] = conj(t[N-i-1])*static_cast<floatingPointType>(i+1);
    }


    for(int i = 0; i < N; i++) {
        w[i].r = s[i].real();
        w[i].i = s[i].imag();
        tx[i].r = t[i].real();
        tx[i].i = t[i].imag();
    }


    for(int i = N; i < 2*N; i++) {
        w[i].r = 0.0;
        tx[i].r = 0.0;
        w[i].i = 0.0;
        tx[i].i = 0.0;
    }


    tx[N].r = t[0].real();
    tx[N].i = t[0].imag();




    kiss_fft_cpx TX[2*N];
    kiss_fft_cpx W[2*N];

    kiss_fft_cpx fft1[2*N];
    kiss_fft(cfg, tx, TX);
    kiss_fft(cfg, w, W);


    for(int i = 0; i < 2*N; i++) {
        fft1[i].r = TX[i].r * W[i].r - TX[i].i * W[i].i;
        fft1[i].i = TX[i].r * W[i].i + TX[i].i * W[i].r;
    }



    t[0] = 0.0;
    for(int i = 1; i < N; i++) {
        t[i] = conj(inputVector[N - i]);
    }
    for(int i = 0; i < N; i++) {
        s[i] = conj(t[N-i-1])*static_cast<floatingPointType>(i+1);
    }
    for(int i = 0; i < N; i++) {
        w[i].r = s[i].real();
        w[i].i = s[i].imag();
        tx[i].r = t[i].real();
        tx[i].i = t[i].imag();
    }
    for(int i = N; i < 2*N; i++) {
        w[i].r = 0.0;
        tx[i].r = 0.0;
        w[i].i = 0.0;
        tx[i].i = 0.0;
    }



    kiss_fft_cpx fft2[2*N];
    kiss_fft(cfg, tx, TX);
    kiss_fft(cfg, w, W);


    tx[N].r = t[0].real();
    tx[N].i = t[0].imag();


    for(int i = 0; i < 2*N; i++) {
        fft2[i].r = TX[i].r * W[i].r - TX[i].i * W[i].i;
        fft2[i].i = TX[i].r * W[i].i + TX[i].i * W[i].r;
    }

    kiss_fft_free(cfg);
    cfg = kiss_fft_alloc(2*N, 1, NULL, NULL);

    kiss_fft_cpx tmp1[2*N];
    kiss_fft_cpx tmp2[2*N];
    kiss_fft(cfg, fft1,tmp1);
    kiss_fft(cfg, fft2,tmp2);


    kiss_fft_free(cfg);




    complex<floatingPointType>* phi = new complex<floatingPointType>[2*N-1];
    for(int i = 0; i < N; i++) {
        phi[i] = complex<floatingPointType>( (tmp1[i].r - tmp2[i].r)/(2*N), ( -tmp1[i].i + tmp2[i].i)/(2*N)    );
    }
    for(int i = 0; i < N-2; i++) {
        phi[N+i] = conj( phi[N-i-2 ] );
    }

    return phi;
}



template <typename floatingPointType>
tuple<complex<floatingPointType>*, floatingPointType, complex<floatingPointType>*> UtilityMathFunctions<floatingPointType>::levinson(const complex<floatingPointType>* inputVector, size_t N) {

    //Taken from https://github.com/cokelaer/spectrum/blob/master/src/spectrum/levinson.py

    int i, j,k,kj;
    int khalf;
    complex<floatingPointType> r[N];
    for(i = 0; i<N; i++) {
        r[i] = inputVector[i];
    }
    floatingPointType T0 = r[0].real();

    complex<floatingPointType> T[N-1];
    for(i = 0; i< N-1; i++) {
        T[i] = r[i+1];
    }
    int M = N - 1;

    complex<floatingPointType>* A = new complex<floatingPointType>[N];
    complex<floatingPointType>* refer = new complex<floatingPointType>[M];
    //fill(begin(A),end(A),0.0);
    //fill(begin(refer),end(refer),0.0);
    A[0] = 1.0;
    for(i = 0; i < M; i++) {
        A[i+1] = 0;
        refer[i] = 0;
    }

    floatingPointType P = T0;
    complex<floatingPointType> save;
    complex<floatingPointType> temp;

    bool warnedSingularMatrix = false;

    for(k = 0; k < M; k++) {
        save = T[k];
        if(k == 0) {
            temp = -save/P;
        } else {
            for(j = 0; j < k; j++) {
                save = save + A[j+1]*T[k-j-1];
            }
            temp = -save/P;
        }
        P = P*(1.0 - (temp.real()*temp.real() + temp.imag()*temp.imag() ));
        if(P < 0 && !warnedSingularMatrix) {
            warnedSingularMatrix = true;
            cout << "Singular matrix " << endl;

            saveArrayToFile(inputVector, N, "D:/data/OffendingVector.txt");


        }
        A[k+1] = temp;
        refer[k] = temp;
        if(k == 0) {
            continue;
        }
        khalf = (k+1)/2;
        for(j = 0; j < khalf; j++) {
            kj = k - j - 1;
            save = A[j+1];
            A[j+1] = save + temp*A[kj+1];
            if(j != kj) {
                A[kj+1] += temp*save;
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

    floatingPointType dx;
    int initialIndex = static_cast<int>(N*x/(splines_[N-1]->x0-splines_[0]->x0));
    if(splines_[initialIndex]->x0 > x) {
        //go left
        for(int i = initialIndex; i>0; i--) {
            dx = fabs(splines_[i]->x0 - splines_[i-1]->x0);
            if( fabs(splines_[i]->x0 - x ) <= dx  ) {
                return splines_[i-1]->evaluate(x);
            }
        }

    } else {
        //go right
        for(int i = initialIndex; i<N-1; i++) {
            dx = fabs(splines_[i]->x0 - splines_[i+1]->x0);
            if( fabs(splines_[i]->x0 - x ) <= dx  ) {
                return splines_[i]->evaluate(x);
            }
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
        splines[i] = new Spline<floatingPointType> {a[i],b[i],c[i],d[i],x[i]};
    }


    SplineInterpolation* interp = new SplineInterpolation(splines,n);


    return  interp;
}



