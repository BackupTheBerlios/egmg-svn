/** \file CycleType.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * 
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
public:
    /**
     * \brief The constructor of a CycleType object
     * 
     * A CylcleType object needs a constructor with at least these parameters:
     * \param[in] maximalDepth  Number of Grid Levels
     */
    CycleType(const Index maximalDepth);
    
    /**
     * \brief The copy constructor of a CycleType object
     * 
     * The copy constructor of a CycleType object has to take care that the
     * copy constructed CycleType object knows everything to do the next
     * coarser grid
     * 
     * \param[in] rhs   the WCycle to copy
     */
    CycleType(const CycleType& rhs);
    
    ~CycleType();
    
    /**
     * \brief solve() on this grid level directly?
     * 
     * solve() checks if we are on the coarsest grid level.
     * 
     * \return  true if we are on the coarsest grid
     */
    bool solve() const;
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    bool repeat();
    
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
    void accelerate(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const Index nx,
        const Index ny);
};

}

#endif /*CYCLETYPE_H_*/
