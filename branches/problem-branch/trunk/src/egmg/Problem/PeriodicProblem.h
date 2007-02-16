#ifndef PERIODICPROBLEM_H_
#define PERIODICPROBLEM_H_

#include "Problem.h"

namespace mg
{

class PeriodicProblem : public mg::Problem
{
public:
	PeriodicProblem(Stencil& stencil, Index nx, Index ny);
	virtual ~PeriodicProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& array ) const;
    virtual DiscreteFunction residuum();
    virtual Point getFirstPoint() const;
    virtual Point getLastPoint() const;
    virtual PeriodicProblem* getCoarsGridProblem(Index nx, Index ny) const;
protected:
    virtual void fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const;
};

}

#endif /*PERIODICPROBLEM_H_*/
