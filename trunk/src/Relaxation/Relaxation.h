/** \file Relaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Relaxation
 */
#ifndef RELAXATION_H_
#define RELAXATION_H_


#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief Relaxation is a 2D Relaxation operator
 * \todo Implement ILU Smoother
 */
class Relaxation
{
private:
	int preSmoothingSteps_;
	int postSmoothingSteps_;
	
public:
	/**
	 * \brief The constructor of a Relaxation object
	 * 
	 * Relaxation constructs a Relaxation object with:
	 * \param[in] preSmoothingSteps		number of pre smoothing steps  (def. 1)
	 * \param[in] postSmoothingSteps	number of post smoothing steps (def. 1)
	 */
	Relaxation(
		const int preSmoothingSteps =1,
		const int postSmoothingSteps =1):
		preSmoothingSteps_(preSmoothingSteps),
		postSmoothingSteps_(postSmoothingSteps){}

    virtual ~Relaxation() {}
    
    /**
	 * \brief getPreSmoothingSteps() returns the number of pre smothing steps
	 * \return the number of pre somthing steps
	 */
	int getPreSmoothingSteps() const {return preSmoothingSteps_;};
	
	/**
	 * \brief setPreSmoothingSteps() sets the number of pre smothing steps
	 * \param[in] preSmoothingSteps		the number of pre somthing steps
	 */
	void setPreSmoothingSteps(int preSmoothingSteps)
	{
		preSmoothingSteps_=preSmoothingSteps;
	}
	
	/**
	 * \brief getPostSmoothingSteps() returns the number of post somthing steps
	 * \return the number of post somthing steps
	 */
	int getPostSmoothingSteps()	const {return postSmoothingSteps_;};
	
	/**
	 * \brief setPostSmoothingSteps() sets the number of post smothing steps
	 * \param[in] postSmoothingSteps		the number of post somthing steps
	 */
	void setPostSmoothingSteps(int postSmoothingSteps)
	{
		postSmoothingSteps_=postSmoothingSteps;
	}
    
    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * \param[in,out] u         the vector representation of the 2D grid to
     *                          perform the relaxation on
     * \param[in] f             the right hand side of the pde
     * \param[in] stencil       the stencil rep. of the pde
     * \param[in] nx            number of steps in x direction
     * \param[in] ny            number of steps in y direction
     */
    virtual void relax(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const=0;
    
    /**
     * \brief does the preSmooth() ing
     * 
     * \param[in,out] u         the vector representation of the 2D grid to
     *                          perform the smoothin on
     * \param[in] f             the right hand side of the pde
     * \param[in] stencil       the stencil rep. of the pde
     * \param[in] nx            number of steps in x direction
     * \param[in] ny            number of steps in y direction
     */    
    virtual void preSmooth(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const
    {
        for(int i=0; i<preSmoothingSteps_; ++i)
            relax(u,f,stencil,nx,ny);
    }
    
    /**
     * \brief does the postSmooth() ing
     * 
     * \param[in,out] u         the vector representation of the 2D grid to
     *                          perform the smoothing on
     * \param[in] f             the right hand side of the pde
     * \param[in] stencil       the stencil rep. of the pde
     * \param[in] nx            number of steps in x direction
     * \param[in] ny            number of steps in y direction
     */   
    virtual void postSmooth(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const
    {
        for(int i=0; i<postSmoothingSteps_; ++i)
            relax(u,f,stencil,nx,ny);
    }
};

}

#endif /*RELAXATION_H_*/
