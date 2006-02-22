/** \file GS_RB.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief GS_RB.h contains the interface of the class GS_RB.
 * \see Relaxation.h
 */
#ifndef GS_RB_H_
#define GS_RB_H_

#include "Relaxation.h"

namespace mg
{

/**
 * \brief GS_RB is a class for a Gauss Seidel relaxation
 * 
 * GS_RB represents a Gauss Seidel relaxation with red-black ordering.
 */
class GS_RB : public mg::Relaxation
{
private:
	int nreld_;
	int nrelu_;
public:
	/**
	 * \brief The constructor of a GS_RB object
	 * 
	 * GS_RB constructs a GS_RB object with:
	 * \param[in] nreld		number of pre smothing steps (default 1)
	 * \param[in] nrelu		number of post smothing steps (default 1)
	 */
	GS_RB(const int nreld =1,const int nrelu =1)
		: nreld_(nreld), nrelu_(nrelu) {}
	virtual ~GS_RB() {}
	/**
	 * \brief getPreSmoothingSteps() returns the number of pre smothing steps
	 * \return the number of pre somthing steps
	 */
	int get_nreld() const {return nreld_;};
	/**
	 * \brief setPreSmoothingSteps() sets the number of pre smothing steps
	 * \param[in] nreld		the number of pre somthing steps
	 */
	void set_nreld(int nreld) {nreld_=nreld;}
	/**
	 * \brief getPostSmoothingSteps() returns the number of post somthing steps
	 * \return the number of post somthing steps
	 */
	int get_nrelu()	const {return nrelu_;};
	/**
	 * \brief setPostSmoothingSteps() sets the number of post smothing steps
	 * \param[in] nrelu		the number of post somthing steps
	 */
	void set_nrelu(int nrelu) {nrelu_=nrelu;}
	/**
	 * \brief relax() executes one relaxation step on the input vector
	 * 
	 * relax() exectues one Gaus seidel relaxation step on the input vector
	 * one a rectangular 2D gird with red-black ordering and the
	 * discretazation given by stencil for a pde.
	 * 
	 * \param[in,out] u		the vector representation of the 2D grid to perform
	 * 						the relaxation on
	 * \param[in] fv		the right hand side of the pde
	 * \param[in] stencil	the stencil rep. of the pde
	 * \param[in] Nx		number of steps in x direction
	 * \param[in] Ny		number of steps in y direction
	 */
	void relax(std::valarray<Precision>& u,
				const std::valarray<Precision>& fv,
				const Stencil& stencil,
				const size_t Nx,const size_t Ny) const;
};

}

#endif /*GS_RB_H_*/
