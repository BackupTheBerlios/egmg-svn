/** \file cycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the template function cycle.
 * \see cycle.h
 */

#include <algorithm>
#include <stdexcept>
#include "cycle.h"
#include "residuum.h"

namespace mg
{
template<typename CycleType>
void cycle(
    CycleType cycleType,
    NumericArray& u,
    const NumericArray& f,
    Stencil& stencil,
    const Prolongation& prolongation,
    const Restriction& restriction,
    Relaxation& relaxation,
    const Index nx,
    const Index ny)
{
    if (cycleType.solve() || std::min(nx,ny)<=6)
    {
        /**
         * \todo implement a direct solver
         */
        for (int i=0; i<10; i++)
            relaxation.relax(u,f,stencil,nx,ny);
    }
    else
    {
        //if it is not possible to do standart coarsening throw an exeption
        if (nx%2!=0 || ny%2!=0)
            throw std::domain_error("u");
        while(cycleType.repeat())
        {
            relaxation.preSmooth(u,f,stencil,nx,ny);
            //calculate the residuum
            NumericArray residv=residuum(u,f,stencil,nx,ny);
            //restrict the residuum to the coars grid
            NumericArray coarsResiduum=restriction.restriction
                    (residv,stencil,prolongation,nx,ny);
            //we going to a coarser grid so Galerkin Operator needs to know
            //the transfer operators
            stencil.pushProlongation(prolongation);
            stencil.pushRestriction(restriction);
            const Index nxNew = nx/2;
            const Index nyNew = ny/2;
            NumericArray coarsGridCor
                    (0.0,(nxNew+1)*(nyNew+1));
            //do a multigrid cycle on the coars grid
            cycle(
                cycleType,
                coarsGridCor,
                coarsResiduum,
                stencil,
                prolongation,
                restriction,
                relaxation,
                nxNew,nyNew);
            //prolongate the coars grid correction to the fine grid
            //approximation
            u+=prolongation.prolongate(coarsGridCor,stencil,nxNew,nyNew);
            //we are going to a smaler grid so remove transfer operators
            stencil.popRestriction();
            stencil.popProlongation();
            relaxation.postSmooth(u,f,stencil,nx,ny);
            cycleType.accelerate(u,f,stencil,nx,ny);
        }
    }
}
}
