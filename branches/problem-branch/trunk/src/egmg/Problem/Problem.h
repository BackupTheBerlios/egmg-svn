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
    virtual const NumericArray& getRightHandSide();
    virtual void setBoundaryConstraint( const Function& boundaryConstraint ) =0;
    virtual void applyBoundaryConstraint() =0;
    virtual void applyBoundaryConstraint( NumericArray& array ) const =0;
    virtual NumericArray& getSolution() =0;
    virtual const NumericArray& getSolution() const =0;
    virtual NumericArray residuum() =0;
    virtual Point getLowerLeftCorner() const =0;
    virtual Problem* getCoarsGridProblem(Index nx, Index ny) const =0;
    virtual Point getFirstPoint() const =0;
    virtual Point getLastPoint() const =0;
    Index getNx() const;
    Index getNy() const;
    
protected:
    Stencil& stencil_;
    const Index nx_;
    const Index ny_;
    NumericArray rightHandSide_;
};

}

#endif /*PROBLEM_H_*/
