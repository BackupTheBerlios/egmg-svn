/** \file CycleType.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the interface CycleType.
 * \see cycle.h
 */
#ifndef CYCLETYPE_H_
#define CYCLETYPE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief CycleType is the interface of a general CycleType
 */
class CycleType
{
private:
    int preSmoothingSteps_;
    int postSmoothingSteps_;

    /**
     * \brief The copy constructor of a CycleType object
     * 
     * The copy constructor of a CycleType object should not be used.
     * \param[in] rhs   the CycleType to copy
     */
    CycleType(const CycleType& rhs);
public:
    /**
     * \brief The constructor of a CycleType object
     * 
     * CycleType constructs a CycleType object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     */
    CycleType(
        const int preSmoothingSteps,
        const int postSmoothingSteps):
        preSmoothingSteps_(preSmoothingSteps),
        postSmoothingSteps_(postSmoothingSteps) {}

    virtual ~CycleType() {}
    
    /**
     * \brief solve() on this grid level directly?
     * 
     * solve() checks if we are on the coarsest grid level.
     * 
     * \return  true if we are on the coarsest grid
     */
    virtual bool solve() const =0;
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    virtual bool repeat() =0;
    
    /**
     * \brief tells a CycleType Object that the next grid level will be entered
     */
    virtual void incrementGridLevel() =0;
    
    /**
     * \brief tells a CycleType Object that a grid level will be left
     */
    virtual void decrementGridLevel() =0;
    
    /**
     * \brief getPreSmoothingSteps() returns the number of pre smothing steps
     * \return the number of pre somthing steps
     */
    int getPreSmoothingSteps() const {return preSmoothingSteps_;};
    
    /**
     * \brief setPreSmoothingSteps() sets the number of pre smothing steps
     * \param[in] preSmoothingSteps     the number of pre somthing steps
     */
    void setPreSmoothingSteps(int preSmoothingSteps)
    {
        preSmoothingSteps_=preSmoothingSteps;
    }
    
    /**
     * \brief getPostSmoothingSteps() returns the number of post somthing steps
     * \return the number of post somthing steps
     */
    int getPostSmoothingSteps() const {return postSmoothingSteps_;};
    
    /**
     * \brief setPostSmoothingSteps() sets the number of post smothing steps
     * \param[in] postSmoothingSteps        the number of post somthing steps
     */
    void setPostSmoothingSteps(int postSmoothingSteps)
    {
        postSmoothingSteps_=postSmoothingSteps;
    }
    
    /**
     * \brief accelerate() the current solution
     * 
     * accelerate can be used to do krylov space accelerations
     * 
     * \param[in,out] u     the vector representation of the 2D grid
     * \param[in] f         the right hand side of the pde
     * \param[in] stencil   the stencil rep. of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     */
    virtual void accelerate(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const Index nx,
        const Index ny) =0;
};

}

#endif /*CYCLETYPE_H_*/
