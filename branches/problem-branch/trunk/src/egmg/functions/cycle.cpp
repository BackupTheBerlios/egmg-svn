/** \file cycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the template function cycle.
 * \see cycle.h
 */

#include <algorithm>
#include <stdexcept>
#include "cycle.h"
#include "residuum.h"
#include "directSolver.h"
#include "printStencil.h"
#include "../general/DiscreteFunction.h"
#include <iostream>
#include <iomanip>

namespace mg
{

void cycle(
    CycleType& cycleType,
    Problem& problem,
    const Prolongation& prolongation,
    const Restriction& restriction,
    const Relaxation& relaxation)
{
    cycleType.incrementGridLevel();
    
    DiscreteFunction& u = problem.getSolution();
    Stencil& stencil = problem.getStencil();
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    
    if ( cycleType.solve())
    {        
        NumericArray temp = directSolver(problem.getLinearEquationSystem());
        problem.setSolution(temp);
    }
    else
    {
        while(cycleType.repeat())
        {
            for(int i=0; i<cycleType.getPreSmoothingSteps(); ++i)
                relaxation.relax(problem);
            //calculate the residuum
            DiscreteFunction residv=problem.residuum();
            //restrict the residuum to the coars grid
            DiscreteFunction coarsResiduum=restriction.restriction
                    (problem,residv);
            const Index nxNew = nx/2;
            const Index nyNew = ny/2;
            //we going to a coarser grid so Galerkin Operator needs to know
            //the transfer operators
            stencil.pushTransferOperators(restriction,prolongation,nxNew,nyNew);
            Problem* coarsGridProblem = problem.getCoarsGridProblem(nxNew,nyNew);
            coarsGridProblem->setRightHandSide( coarsResiduum );
            //do a multigrid cycle on the coars grid

            cycle(
                cycleType,
                *coarsGridProblem,
                prolongation,
                restriction,
                relaxation);

            //prolongate the coars grid correction to the fine grid
            //approximation
            u+=prolongation.prolongate(*coarsGridProblem);
            
            //we are going to a smaler grid so remove transfer operators
            stencil.popTransferOperators();

            for(int i=0; i<cycleType.getPostSmoothingSteps(); ++i)
                relaxation.relax(problem);

            cycleType.accelerate(problem);
            
            delete coarsGridProblem;
        }
    }

    cycleType.decrementGridLevel();
}
}
