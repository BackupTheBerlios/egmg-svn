#ifndef PERIODICPROBLEM_H_
#define PERIODICPROBLEM_H_

#include "Problem.h"

namespace mg
{

class PeriodicProblem : public mg::Problem
{
public:
	PeriodicProblem(
        Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy);
	virtual ~PeriodicProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& array ) const;
    virtual DiscreteFunction residuum();
    virtual IndexPair getFirstPoint(Index nx, Index ny) const;
    virtual IndexPair getLastPoint(Index nx, Index ny) const;
    virtual ProblemPtr getCoarsGridProblem(
        Index nxNew,
        Index nyNew,
        Precision hxNew,
        Precision hyNew) const;
protected:
    virtual void fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const;
};

}

#endif /*PERIODICPROBLEM_H_*/
