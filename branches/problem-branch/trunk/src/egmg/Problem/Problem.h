#ifndef PROBLEM_H_
#define PROBLEM_H_


#include "../general/parameters.h"
#include "../Function/Function.h"
#include "../Stencil/Stencil.h"
#include "../general/DiscreteFunction.h"

namespace mg
{

class Problem
{
public:
    Problem(Stencil& stencil, Index nx, Index ny);
	virtual ~Problem();
    Stencil& getStencil();
    void setRightHandSide( const Function& rightHandSide );
    void setRightHandSide( const DiscreteFunction& rightHandSide );
    const DiscreteFunction& getRightHandSide();
    virtual void setBoundaryConstraint( const Function& boundaryConstraint ) =0;
    virtual void applyBoundaryConstraint() =0;
    virtual void applyBoundaryConstraint( DiscreteFunction& function ) const =0;
    DiscreteFunction& getSolution();
    const DiscreteFunction& getSolution() const;
    virtual DiscreteFunction residuum() =0;
    virtual Problem* getCoarsGridProblem(Index nx, Index ny) const =0;
    virtual Point getFirstPoint() const =0;
    virtual Point getLastPoint() const =0;
    Index getNx() const;
    Index getNy() const;
    
protected:
    Stencil& stencil_;
    const Index nx_;
    const Index ny_;
    DiscreteFunction rightHandSide_;
    DiscreteFunction solution_;
};

}

#endif /*PROBLEM_H_*/
