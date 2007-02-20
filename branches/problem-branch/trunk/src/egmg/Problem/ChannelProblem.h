#ifndef CHANNELPROBLEM_H_
#define CHANNELPROBLEM_H_

#include "Problem.h"

namespace mg
{

class ChannelProblem : public mg::Problem
{
public:
	ChannelProblem(Stencil& stencil, Index nx, Index ny);
	virtual ~ChannelProblem();
    
    virtual void setBoundaryConstraint( const Function& boundaryConstraint );
    virtual void applyBoundaryConstraint();
    virtual void applyBoundaryConstraint( DiscreteFunction& array ) const;
    virtual DiscreteFunction residuum();
    virtual Point getFirstPoint() const;
    virtual Point getLastPoint() const;
    virtual ChannelProblem* getCoarsGridProblem(Index nx, Index ny) const;
protected:
    virtual void fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const;
};

}

#endif /*CHANNELPROBLEM_H_*/
