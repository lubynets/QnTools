#ifndef QNTOOLS_H1F_H_
#define QNTOOLS_H1F_H_

#include "TH1F.h"

namespace Qn {
  
class H1F : public TH1F {
 public:
  friend H1F Merge(const H1F &, const H1F &);
  
    /// \cond CLASSIMP
 ClassDef(H1F, 1);
  /// \endcond
};
H1F Merge(const H1F &lhs, const H1F &rhs);
}
#endif //QNTOOLS_H1F_H_