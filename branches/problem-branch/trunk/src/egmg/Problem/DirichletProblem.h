#ifndef DIRICHLETPROBLEM_H_
#define DIRICHLETPROBLEM_H_

#include "Problem.h"

namespace mg
{

class DirichletProblem : public mg::Problem
{
public:
	DirichletProblem(
        Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy);
	virtual ~DirichletProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& function ) const;
    virtual DiscreteFunction residuum();
    virtual IndexPair getFirstPoint(Index nx, Index ny) const;
    virtual IndexPair getLastPoint(Index nx, Index ny) const;
    virtual ProblemPtr getCoarsGridProblem(
        Index nxNew, Index nyNew,
        Precision hxNew, Precision hyNew) const;
protected:
    virtual void fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const;
};

}

#endif /*DIRICHLETPROBLEM_H_*/
