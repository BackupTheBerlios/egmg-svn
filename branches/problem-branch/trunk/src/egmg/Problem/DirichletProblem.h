#ifndef DIRICHLETPROBLEM_H_
#define DIRICHLETPROBLEM_H_

#include "Problem.h"

namespace mg
{

class DirichletProblem : public mg::Problem
{
public:
	DirichletProblem(Stencil& stencil, Index nx, Index ny);
	virtual ~DirichletProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& function ) const;
    virtual DiscreteFunction residuum();
    virtual Point getFirstPoint() const;
    virtual Point getLastPoint() const;
    virtual DirichletProblem* getCoarsGridProblem(Index nx, Index ny) const;
};

}

#endif /*DIRICHLETPROBLEM_H_*/
