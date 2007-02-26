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
    Problem(
        Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy);
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
    virtual ProblemPtr getCoarsGridProblem(
        Index nxNew, Index nyNew,
        Precision hxNew, Precision hyNew) const =0;
    IndexPair getFirstPoint() const;
    IndexPair getLastPoint() const;
    virtual IndexPair getFirstPoint(Index nx, Index ny) const =0;
    virtual IndexPair getLastPoint(Index nx, Index ny) const =0;
    Index getNx() const;
    Index getNy() const;
    Precision getHx() const;
    Precision getHy() const;
    Point getOrigin() const;
    LinearEquationSystem getLinearEquationSystem() const;
    void setSolution(NumericArray& solution);
    
protected:
    virtual void fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const =0;
    Stencil& stencil_;
    DiscreteFunction rightHandSide_;
    DiscreteFunction solution_;
    
private:
    void fillRightHandSide(NumericArray& rightSide) const;
    void pointFillMatrix(
        NumericArray& matrix,
        const Index dimension,
        const Position postion,
        const Index sx,
        const Index sy) const;
    void fillMatrix(
        NumericArray& matrix,
        Index dimension) const;
};

}

#endif /*PROBLEM_H_*/
