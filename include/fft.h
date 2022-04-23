
#ifndef FFT
#define FFT

#include <iostream>
#include <cmath>
#include <new>
#include <vector>
#include <complex>
#include "TMath.h"
// #include <numbers>

namespace fft {

  class dft {

  public:

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> CNArray;

    dft(const CNArray& x) {
      arraySize = x.size();
      outputArray.reserve(arraySize);
      outputArray.assign(x.begin(), x.end());
    };
    virtual ~dft() {};

    CNArray getDFT() {copyHalfArray(); return halfOutputArray;};
    CNArray getFullDFT() {return outputArray;};
    
  private:
  
    CNArray outputArray;
    CNArray halfOutputArray;
    int arraySize;
    
    void  DFT(CNArray& x);
    void iDFT(CNArray& x);

    void  doDFT() {DFT(outputArray);};
    void doiDFT() {iDFT(outputArray);};
    
    void copyHalfArray();

  };


}


#endif
