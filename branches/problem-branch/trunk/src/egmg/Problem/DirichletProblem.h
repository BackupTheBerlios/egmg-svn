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
    virtual void applyBoundaryConstraint( NumericArray& array ) const;
    virtual NumericArray residuum();
    virtual NumericArray& getSolution();
    virtual const NumericArray& getSolution() const;
    virtual Point getLowerLeftCorner() const;
    virtual Point getUpperRightCorner() const;
    virtual Point getFirstPoint() const;
    virtual Point getLastPoint() const;
    virtual DirichletProblem* getCoarsGridProblem(Index nx, Index ny) const;
    
private:
    NumericArray solution_;
};

}

#endif /*DIRICHLETPROBLEM_H_*/
