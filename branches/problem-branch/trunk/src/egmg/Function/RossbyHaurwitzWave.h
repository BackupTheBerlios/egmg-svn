#ifndef ROSSBYHAURWITZWAVE_H_
#define ROSSBYHAURWITZWAVE_H_


#include "Function.h"

namespace mg
{

class RossbyHaurwitzWave: public Function
{
public:
    RossbyHaurwitzWave(
        Precision t=0.0,
        Precision c=1.0,
        Precision beta=1.0,
        Precision k=2.0*Pi,
        Precision l=2.0*Pi);
    void setT(Precision t);
    Precision getT();
    Precision operator() (Precision x, Precision y, Precision t) const;
protected:
	virtual Precision apply(Precision x, Precision y) const;
private:
    Precision t_;
    const Precision c_;
    const Precision beta_;
    const Precision k_;
    const Precision l_;
};

}

#endif /*ROSSBYHAURWITZWAVE_H_*/
