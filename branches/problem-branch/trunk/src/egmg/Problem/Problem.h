#ifndef PROBLEM_H_
#define PROBLEM_H_


#include "../general/parameters.h"
#include "../Function/Function.h"
#include "../Stencil/Stencil.h"

namespace mg
{

class Problem
{
public:
    Problem(Stencil& stencil, Index nx, Index ny);
	virtual ~Problem();
    virtual Stencil& getStencil();
    virtual void setRightHandSide( const Function& rightHandSide );
    virtual void setRightHandSide( const NumericArray& rightHandSide );
    virtual const NumericArray& NumericArray getRightHandSide();
    virtual void setBoundaryConstraint( const Function& boundaryConstraint ) =0;
    virtual void applyBoundaryConstraint() =0;
    virtual NumericArray& getSolution() =0;
    virtual NumericArray residuum() =0;
    
protected:
    const Index nx_;
    const Index ny_;
private:
    Stencil& stencil_;
    NumericArray rightHandSide_;
};

}

#endif /*PROBLEM_H_*/
