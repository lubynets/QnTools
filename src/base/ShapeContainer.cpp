#include "ShapeContainer.hpp"

// ClassImp(ShapeContainer);    // linking error occurs

namespace Qn {

void ShapeContainer::SetShape(TH1F* histosgnl, TH1F* histobckgr)
{
  histo_sgnl_ = histosgnl;
  histo_bckgr_ = histobckgr;
  is_histo_sgnl_ = true;
  is_histo_bckgr_ = true;
}

void ShapeContainer::SetShape(TF1* funcsgnl, TF1* funcbckgr)
{
  func_sgnl_ =  funcsgnl;
  func_bckgr_ = funcbckgr;
  is_histo_sgnl_ = false;
  is_histo_bckgr_ = false;
}

void ShapeContainer::SetShape(TH1F* histosgnl, TF1* funcbckgr)
{
  histo_sgnl_ = histosgnl;
  func_bckgr_ = funcbckgr;
  is_histo_sgnl_ = true;
  is_histo_bckgr_ = false;
}

void ShapeContainer::SetShape(TF1* funcsgnl, TH1F* histobckgr)
{
  func_sgnl_ =  funcsgnl;
  histo_bckgr_ = histobckgr;
  is_histo_sgnl_ = false;
  is_histo_bckgr_ = true;
}

float ShapeContainer::GetSignal(float x)                       //TODO add warnings or even errors if histo or/and function is not initialized
{
  if(is_histo_sgnl_ == true)
    return histo_sgnl_ -> Interpolate(x);
  else
    return func_sgnl_ -> Eval(x);
}

float ShapeContainer::GetBackground(float x)                   //TODO add warnings or even errors if histo or/and function is not initialized
{
  if(is_histo_bckgr_ == true)
    return histo_bckgr_ -> Interpolate(x);
  else
    return func_bckgr_ -> Eval(x);
}

float ShapeContainer::GetSignalIntegral(float left, float right)
{
  if(is_histo_sgnl_ == true)
    return ShapeContainer::HistoIntegral(histo_sgnl_, left, right);
  else
    return func_sgnl_ -> Integral(left, right) / bin_width_;
}

float ShapeContainer::GetBackgroundIntegral(float left, float right)
{
  if(is_histo_bckgr_ == true)
    return ShapeContainer::HistoIntegral(histo_bckgr_, left, right);
  else
    return func_bckgr_ -> Integral(left, right) / bin_width_;
}

float ShapeContainer::HistoIntegral(TH1F* histo, float low, float up)
{
  if(up<low)
    return -ShapeContainer::HistoIntegral(histo, up, low);
  
  if(low<histo->GetXaxis()->GetBinLowEdge(0))
    low = histo->GetXaxis()->GetBinLowEdge(0);
  
  if(up>histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1))
    up = histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1);
  
  const int lowbin = histo->GetXaxis()->FindBin(low);
  const int upbin = histo->GetXaxis()->FindBin(up);
  
  float integral = histo->Integral(lowbin, upbin);
  integral -= histo->GetBinContent(lowbin)*(low-histo->GetXaxis()->GetBinLowEdge(lowbin))/histo->GetXaxis()->GetBinWidth(lowbin);
  integral -= histo->GetBinContent(upbin)*(histo->GetXaxis()->GetBinLowEdge(upbin+1)-up)/histo->GetXaxis()->GetBinWidth(upbin);
  
  return integral;
}

// float ShapeContainer::FuncIntegral(TF1* func, float low, float up) {
//   const int nbins = 1000;
//   const float binwidth = (up-low)/nbins;
//   TH1F* histo = new TH1F("h1", "", nbins, low, up);
//   for(int i=1; i<=nbins; i++)
//     histo->SetBinContent(i, func->Eval(low+(i-0.5)*binwidth));
//   
//   return HistoIntegral(histo, low, up)*binwidth;
// }

}