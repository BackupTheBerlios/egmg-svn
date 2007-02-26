#ifndef CHANNELPROBLEM_H_
#define CHANNELPROBLEM_H_

#include "Problem.h"

namespace mg
{

class ChannelProblem : public mg::Problem
{
public:
	ChannelProblem(
        Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy);
	virtual ~ChannelProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& array ) const;
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

#endif /*CHANNELPROBLEM_H_*/
