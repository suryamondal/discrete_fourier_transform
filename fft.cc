
#include "fft.h"

void fft::dft::FFT(fft::CNArray& x) {
  
  const int N = x.size();
  if (N <= 1) return;

  // divide
  fft::CNArray even, odd;
  even.clear(); odd.clear();
  for (Long64_t k = 0; k < N; k += 2) {
    if(k+1 == N) {break;}
    even.push_back(x[k]);
    odd.push_back(x[k+1]);
  }
  
  // compute
  fft::dft::doFFT(even);
  fft::dft::doFFT(odd);

  // combine
  for (Long64_t k = 0; k < N/2; k++) {
    Complex t = std::polar(1.0, - 2. * Pi() * k / N) * odd[k];
    x[k    ] = even[k] + t;
    x[k+N/2] = even[k] - t;
  }  
}

void fft::dft::iFFT(fft::CNArray& x) {

  const int N = x.size();

  // conjugate the complex numbers
  for (Long64_t k=0;k<N;k++) {
    x[k] = std::conj(x[k]);
  }

  // forward fft
  fft::dft::FFT(x);

  // conjugate the complex numbers again
  for (Long64_t k=0;k<N;k++) {
    x[k] = std::conj(x[k]);
  }

  // scaling
  for (Long64_t k=0;k<N;k++) {
    x[k] /= N;
  }
}

void fft::dft::doFFT() {
  fft::dft::FFT(fft::dft::outputArray);
}

void fft::dft::doiFFT() {
  fft::dft::iFFT(fft::dft::outputArray);
}

void fft::dft::copyHalfArray() {
  fft::dft::halfOutputArray.reserve(fft::dft::arraySize/2);
  fft::dft::halfOutputArray.assign(outputArray.begin(),
				   (outputArray.begin() +
				    fft::dft::arraySize/2 - 1));
}

fft::CNArray fft::dft::GetFFT() {
  fft::dft::copyHalfArray();
  return fft::dft::halfOutputArray;
}

