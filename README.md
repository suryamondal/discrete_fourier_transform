# Discrete Fourier Transform

This is a basic code for fast Fourier transform based on Cooley-Turkey.
Here, a few a combination of frequency is generated. Then Fourier
transform is performed on the spectrum.

A few frequencies then removed from the Fourier spectra. The inverse
Fourier transform is performed on the Fourier spectra.

Compile and run with:
```
mkdir -p build
cd build
make
cd ..
./fft_main
```

For more info, one might check this [website](https://fakephysicist.com/misc/cooley-turkey-dfft-vectors-cpp/).