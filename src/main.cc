
#include "TPad.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"

#include "utils.h"
#include "fft.h"


const double tPD = 20.;		 // Total Path Difference in sec
const double maxFrequency = 10.; // Maximum Frequency to Calculate

const Long64_t maxGenPt = 1.;	     // If you force
const Long64_t maxDataPt = std::pow(2,1); // If you force

int main() {

  const Long64_t minPt = 3.*(maxFrequency*tPD+1);
  const Long64_t inputDataPt = maxGenPt*minPt;

  std::cout << " Generating Data...." << std::endl;

  std::vector<double> pos0,amp0;	// Input Data Vector before Padding
  pos0.clear(); pos0.reserve(inputDataPt);
  amp0.clear(); amp0.reserve(inputDataPt);
  
  /** Generating Data Points */
  for(int ij=0;ij<inputDataPt;ij++) {
    double tt = ij*tPD/inputDataPt;
    double amp =
      1.*sin(2.*TMath::Pi()*4.5*tt) +
      2.*sin(2.*TMath::Pi()*4.6*tt) +
      3.*sin(2.*TMath::Pi()*4.7*tt) +
      4.*sin(2.*TMath::Pi()*4.8*tt) +
      5.*sin(2.*TMath::Pi()*4.9*tt) +
      6.*sin(2.*TMath::Pi()*5.0*tt) +
      7.*sin(2.*TMath::Pi()*5.1*tt) +
      8.*sin(2.*TMath::Pi()*5.2*tt) +
      9.*sin(2.*TMath::Pi()*5.3*tt) +
      10.*sin(2.*TMath::Pi()*5.4*tt) ;
    pos0.push_back(tt);
    amp0.push_back(amp);
  }
  
  std::cout << " Input Data Size " << inputDataPt << std::endl;

  const Long64_t maxDesiredPt = std::max(inputDataPt,maxDataPt);

  int pFactor = utils::findExponent(maxDesiredPt,2);
  const Long64_t nn = std::pow(2,pFactor); // Number for Padding
  std::cout << " Size of FFT Input " << nn << std::endl;

  std::vector<double> outPos,outAmp; // Vector for Padded Data
  outPos.clear(); outPos.reserve(nn);
  outAmp.clear(); outAmp.reserve(nn);

  utils::padData(pos0, amp0, outPos, outAmp, nn);

  fft::dft::CNArray data_C;		// vector for Complex number
  data_C.clear(); data_C.reserve(nn);
  for(int ij=0;ij<nn;ij++) {
    fft::dft::Complex tempComp(outAmp[ij],0.); // Setting the Real Part
    data_C.push_back(tempComp);
  }

  /** d-FFT */

  std::cout << " DFT Starts...." << std::endl;
  fft::dft fft_straight(data_C);
  fft_straight.doDFT();
  fft::dft::CNArray out_data_C = fft_straight.getDFT();
  std::cout << " DFT Ends...." << std::endl;

  double maxHeight = 0;
  std::vector<double> cBin,cAmp;
  cBin.clear(); cBin.reserve(nn/2);
  cAmp.clear(); cAmp.reserve(nn/2);
  for(int ij=0;ij<nn/2;ij++) {
    cBin.push_back(ij/tPD);
    double tmpAmp = abs(out_data_C[ij]);
    cAmp.push_back(tmpAmp);
    if(maxHeight<tmpAmp) {maxHeight = tmpAmp;}
  }

  /** Inverse d-FFT */

  fft::dft::CNArray data_I = fft_straight.getFullDFT();

  /** removing background bellow 50% */
  for(auto it=data_I.begin(); it!=data_I.end(); it++) {
    if(abs(*it)<maxHeight*0.5) {
      fft::dft::Complex tempComp(0.,0.);
      *it = tempComp; }}
  
  std::cout << " iDFFT Starts...." << std::endl;
  fft::dft fft_inverse(data_I);
  fft_inverse.doiDFT();
  fft::dft::CNArray out_data_I = fft_inverse.getFullDFT();
  std::cout << " iDFFT Ends...." << std::endl;

  std::vector<double> iBin,iAmp;
  iBin.clear();iAmp.clear();
  for(int ij=0;ij<nn;ij++) {
    iBin.push_back(ij*tPD/(nn-1.));
    iAmp.push_back(real(out_data_I[ij]));
  }
    
  TCanvas *c1 = new TCanvas("c1","c1",2400,2400);
  c1->Divide(1,4);

  TGraph* h3 = new TGraph(inputDataPt,&pos0[0],&amp0[0]); // Gen
  TGraph* h0 = new TGraph(nn,&outPos[0],&outAmp[0]);	  // Pad
  TGraph* h1 = new TGraph(int(nn/2),&cBin[0],&cAmp[0]);	  // d-FFT
  TGraph* h2 = new TGraph(nn,&iBin[0],&iAmp[0]);	  // id-FFT
  
  c1->cd(1);
  h3->SetTitle("Generated");
  h3->Draw("AL*");
  h3->GetXaxis()->SetRangeUser(0., tPD);
  c1->cd(2);
  h0->SetTitle("Padded");
  h0->Draw("AL");
  h0->GetXaxis()->SetRangeUser(0., tPD);
  c1->cd(3);
  h1->SetTitle("DFT");
  h1->Draw("AL*");
  // gPad->SetLogy(1);
  h1->GetXaxis()->SetRangeUser(4., 6.);
  c1->cd(4);
  h2->SetTitle("Inverse transform");
  h2->Draw("AL");
  h2->GetXaxis()->SetRangeUser(0., tPD);

  c1->SaveAs("fft.pdf");
  // c1->SaveAs("fft.C");
  
  return 0;
}
