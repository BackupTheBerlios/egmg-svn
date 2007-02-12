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
    virtual NumericArray residuum();
    
private:
    NumericArray solution_;
};

}

#endif /*DIRICHLETPROBLEM_H_*/
