/** \file GS_Lex.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief GS_Lex.h contains the interface of the class GS_Lex.
 * \see Relaxation.h
 */
#ifndef GS_LEX_H_
#define GS_LEX_H_

#include "Relaxation.h"

namespace mg
{

/**
 * \brief GS_Lex is a class template for a Gauss Seidel relaxation
 * 
 * GS_Lex represents a Gauss Seidel relaxation with the discretazation given by
 * Stencil and lexicographic ordering. Stencil is a template parameter and not a
 * normal class member because of prefomance issues. stencil.apply() is the most
 * frequently called function. To use stencil.apply() as inline functions it is
 * given as a template parameter.
 */
class GS_Lex : public mg::Relaxation
{
private:
	int nreld_;
	int nrelu_;
public:
	/**
	 * \brief The constructor of a GS_Lex object
	 * 
	 * GS_Lex constructs a GS_Lex object with:
	 * \param[in] nreld		number of pre smothing steps (default 1)
	 * \param[in] nrelu		number of post smothing steps (default 1)
	 */
	GS_Lex(const int nreld =1,const int nrelu =1)
		: nreld_(nreld), nrelu_(nrelu) {}
	virtual ~GS_Lex() {}
	
	/**
	 * \brief getPreSmoothingSteps() returns the number of pre smothing steps
	 * \return the number of pre somthing steps
	 */
	int get_nreld() const {return nreld_;}
	
	/**
	 * \brief setPreSmoothingSteps() sets the number of pre smothing steps
	 * \param[in] nreld		the number of pre somthing steps
	 */
	void set_nreld(int nreld) {nreld_=nreld;}
	
	/**
	 * \brief getPostSmoothingSteps() returns the number of post somthing steps
	 * \return the number of post somthing steps
	 */
	int get_nrelu() const {return nrelu_;}
	
	/**
	 * \brief setPostSmoothingSteps() sets the number of post smothing steps
	 * \param[in] nrelu		the number of post somthing steps
	 */
	void set_nrelu(int nrelu) {nrelu_=nrelu;}
	
	/**
	 * \brief relax() executes one relaxation step on the input vector
	 * 
	 * relax() exectues one Gaus seidel relaxation step on the input vector
	 * one a rectangular 2D gird with lexicographic ordering and the
	 * discretazation given by stencil for a pde.
	 * 
	 * \param[in,out] u			the vector representation of the 2D grid to
	 * 							perform the relaxation on
	 * \param[in] fv			the right hand side of the pde
	 * \param[in] stencil		the stencil rep. of the pde
	 * \param[in] Nx			number of steps in x direction
	 * \param[in] Ny			number of steps in y direction
	 */
	void relax(std::valarray<Precision>& u, const std::valarray<Precision>& fv,
				const Stencil& stencil,
				const size_t Nx,const size_t Ny) const;
};

}
#endif /*GS_LEX_H_*/
