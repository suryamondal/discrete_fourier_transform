
#ifndef FFT
#define FFT

#include <iostream>
#include <cmath>
#include <new>
#include <vector>
#include <complex>

namespace fft {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> CNArray;
  
  class dft {

  private:

    CNArray outputArray;
    CNArray halfOutputArray;
    int arraySize;

    dft(const CNArray &x) {
      arraySize = x.size();
      outputArray.reserve(arraySize);
      outputArray.assign(x.begin(), x.end());
    };
    ~dft();

    void  FFT(CNArray &x);
    void iFFT(CNArray &x);

    void  doFFT();
    void doiFFT();

    void copyHalfArray();
    CNArray getFFT();

  };


}


#endif
