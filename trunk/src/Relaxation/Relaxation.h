/** \file Relaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Relaxation
 */
#ifndef RELAXATION_H_
#define RELAXATION_H_
#include<valarray>
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief Relaxation is a 2D Relaxation operator
 * \todo 	move implementaion of get_nreld, set_nreld, get_nrelu and set_nrelu
 * 			into Relaxation
 */
class Relaxation
{
public:
	virtual ~Relaxation() {}
	/**
	 * \brief get_nreld() returns the number of pre smothing steps
	 * \return the number of pre somthing steps
	 */
	virtual int get_nreld() const =0;
	/**
	 * \brief set_nreld() sets the number of pre smothing steps
	 * \param[in] nreld		the number of pre somthing steps
	 */
	virtual void set_nreld(int nreld) =0;
	/**
	 * \brief get_nrelu() returns the number of post somthing steps
	 * \return the number of post somthing steps
	 */
	virtual int get_nrelu()	const =0;
	/**
	 * \brief set_nrelu() sets the number of post smothing steps
	 * \param[in] nrelu		the number of post somthing steps
	 */
	virtual void set_nrelu(int nrelu) =0;
	
	/**
	 * \brief relax() executes one relaxation step on the input vector
	 * 
	 * \param[in,out] u			the vector representation of the 2D grid to
	 * 							perform the relaxation on
	 * \param[in] fv			the right hand side of the pde
	 * \param[in] stencil		the stencil rep. of the pde
	 * \param[in] Nx			number of steps in x direction
	 * \param[in] Ny			number of steps in y direction
	 */
	virtual void relax(std::valarray<precision>& u,
				const std::valarray<precision>& fv,
				const Stencil& stencil,
				const size_t Nx, const size_t Ny) const=0;
};

}

#endif /*RELAXATION_H_*/
