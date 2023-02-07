#ifndef ShapeContainer_H
#define ShapeContainer_H

#include "TH1.h"
#include "TF1.h"
#include "TObject.h"

namespace Qn {

class ShapeContainer : public TObject
{
public:
  
  ShapeContainer() = default;
  virtual ~ShapeContainer() = default;

  void SetInputHisto(TH1F* histo);
  
  void SetShape(TH1F* histosgnl, TH1F* histobckgr);
  void SetShape(TF1* funcsgnl, TF1* funcbckgr);
  void SetShape(TH1F* histosgnl, TF1* funcbckgr);
  void SetShape(TF1* funcsgnl, TH1F* histobckgr);
//   void SetChi2BckgrFit(float value) { chi2_bckgr_fit_ = value ; };
  void SetBinWidth(float value) { bin_width_ = value; };

  float GetSignal(float x);
  float GetBackground(float x);
//   float GetChi2BckgrFit() { return chi2_bckgr_fit_ ; };
  
  TH1F* GetInputHisto() const { return histo_in_; };

  TH1F* GetHistoSgnl() const { return histo_sgnl_; };
  TH1F* GetHistoBckgr() const { return histo_bckgr_; };
  TF1* GetFuncSgnl() const { return func_sgnl_; };
  TF1* GetFuncBckgr() const { return func_bckgr_; };
  
  float GetSignalIntegral(float left, float right);
  float GetBackgroundIntegral(float left, float right);
  
private:

  TH1F* histo_in_{nullptr};
  
  TH1F* histo_sgnl_{nullptr};
  TH1F* histo_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  double bin_width_{0};
//   float chi2_bckgr_fit_{-999.};
  
  bool is_histo_sgnl_{false};
  bool is_histo_bckgr_{false};
  
  float HistoIntegral(TH1F* histo, float low, float up);
//   float FuncIntegral(TF1* histo, float low, float up);
  
  ClassDef(ShapeContainer, 1);
};
}
#endif//ShapeContainer_H
