/** \file cycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>, Matthias Rettenmeier
 * \brief Contains the implimentation of the function cycle.
 * \see cycle.h
 */

#include <algorithm>
#include <stdexcept>
#include "cycle.h"
#include "residuum.h"

namespace mg
{
    void cycle(
        std::valarray<Precision>& u,
        const std::valarray<Precision>& f,
        Stencil& stencil,
        const Prolongation& prolongation,
        const Restriction& restriction,
        Relaxation& relaxation,
        const size_t nx,
        const size_t ny,
        size_t gamma,
        size_t l,
        size_t mode)
    {
        // levelgamma standard set to gamma;
        size_t levelgamma=gamma;

        // if F-cycle ajustments have to be made
        if (mode==0)
        {
            levelgamma=1;
            if (gamma==0)
            {
                mode=2;
                gamma=2;
            }
            else
                mode=1;
        }
        if (std::min(nx,ny)>std::pow(2.0,static_cast<int>(l)))
        {
            //if it is not possible to do standart coarsening throw an exeption
            if (nx%2!=0 || ny%2!=0)
                throw std::domain_error("u");
            for (size_t j=1; j<=levelgamma; j++)
            {
                relaxation.preSmooth(u,f,stencil,nx,ny);
                //calculate the residuum
                std::valarray<Precision> residv=residuum(u,f,stencil,nx,ny);
                //restrict the residuum to the coars grid
                std::valarray<Precision> coarsResiduum=restriction.restriction
                        (residv,stencil,prolongation,nx,ny);
                //we going to a coarser grid so Galerkin Operator needs to know
                //the transfer operators
                stencil.pushProlongation(prolongation);
                stencil.pushRestriction(restriction);
                const size_t nxNew = nx/2;
                const size_t nyNew = ny/2;
                std::valarray<Precision> coarsGridCor
                        (0.0,(nxNew+1)*(nyNew+1));
                //do a multigrid cycle on the coars grid
                cycle(coarsGridCor,coarsResiduum,stencil,prolongation,
                        restriction,relaxation,nxNew,nyNew,l,gamma,mode);
                //if mode one selected call cylce with gamma=1 in next iteration
                if (mode==2)
                    gamma=1;
                //prolongate the coars grid correction to the fine grid
                //approximation
                u+=prolongation.prolongate(coarsGridCor,stencil,nxNew,nyNew);
                //we are going to a smaler grid so remove transfer operators
                stencil.popRestriction();
                stencil.popProlongation();
                relaxation.postSmooth(u,f,stencil,nx,ny);
            }
        }
        else
        {
            /**
             * \todo implement a direct solver
             */
            for (int i=0; i<std::pow(2.0,static_cast<int>(l)); i++)
                relaxation.relax(u,f,stencil,nx,ny);
        }
    }
}
