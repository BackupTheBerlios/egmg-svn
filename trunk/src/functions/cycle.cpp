/** \file cycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>, Matthias Rettenmeier
 * \brief Contains the implimentation of the function cycle.
 * \see cycle.h
 */
#include "cycle.h"
#include "residuum.h"
#include<algorithm>
#include<stdexcept>

namespace mg
{
	void cycle(std::valarray<Precision>& u,const std::valarray<Precision>& fv,
				Stencil& stencil,
				const Prolongation& prolong, const Restriction& restriction,
				Relaxation& relax, const size_t Nx, const size_t Ny,
	    		size_t gamma,size_t l,size_t mode)
	{
		// levelgamma standard set to gamma;
		size_t levelgamma = gamma;

		// if F-cycle ajustments have to be made
		if (mode == 0)
		{
			levelgamma = 1;
			if (gamma == 0)
			{
				mode=2;
          		gamma = 2;
          	}
        	else
        	{
          		mode=1;
        	}
      	}
		if ( std::min(Nx,Ny) > std::pow(2.0,static_cast<int>(l)))
		{
        	//if it is not possible to do standart coarsening throw an exeption
			if ((Nx%2!=0) || ((Ny%2)!=0))
				throw std::domain_error("u");
			for(size_t j=1; j<= levelgamma; j++)
			{
          		//presomthing
				for(int i=0;i<relax.getPreSmoothingSteps();i++)
					relax.relax(u,fv,stencil,Nx,Ny);
				//calculate the residuum
				std::valarray<Precision> residv = residuum(u,fv,stencil,Nx,Ny);
				//restrict the residuum to the coars grid
				std::valarray<Precision> coars_residuum
						= restriction.restriction(residv,stencil,prolong,Nx,Ny);
				//we going to a coarser grid so Galerkin Operator needs to know
				//the transfer operators
				stencil.pushProlongation(prolong);
				stencil.pushRestriction(restriction);
				const size_t Nx_new = Nx/2;
				const size_t Ny_new = Ny/2;
				std::valarray<Precision> coars_grid_cor(
													0.0,(Nx_new+1)*(Ny_new+1));
          		//do a multigrid cycle on the coars grid
				cycle(coars_grid_cor,coars_residuum,stencil,prolong,
						restriction,relax,Nx_new,Ny_new,l,gamma, mode);
				//if mode one selected call cylce with gamma=1 in next iteration
          		if (mode == 2 )
					gamma =1;
				//prolongate the coars grid correction to the fine grid
				//approximation
				u+= prolong.prolongate(coars_grid_cor,stencil,Nx_new,Ny_new);
				//we are going to a smaler grid so remove transfer operators
				stencil.popRestriction();
				stencil.popProlongation();
				//postsomthing
				for(int i=0;i<relax.getPostSmoothingSteps();i++)
					relax.relax(u,fv,stencil,Nx,Ny);
			}
		}
		else
      	{
      		/**
			 * \todo implement a direct solver
			 */
			for (int i=0;i<std::pow(2.0,static_cast<int>(l));i++)
				relax.relax(u,fv,stencil,Nx,Ny);
		}
	}
}
