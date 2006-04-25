/** \file Injection.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the class Injection
 */
#ifndef INJECTION_H_
#define INJECTION_H_

#include "Restriction.h"

namespace mg
{

/**
 * \brief Injection is a 2D restriction operator
 * 
 * Injection is a representation of a 2D restriction operator that uses
 * injection to do its job.
 */
class Injection : public mg::Restriction
{
private:
    const Precision weight_;
    const NumericArray i_;
    const PositionArray jx_;
    const PositionArray jy_;
    //initilize i_, makes it possible to make i_ const
    NumericArray initI_(Precision weight) const
    {
        return NumericArray(weight,1);
    }
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJX_() const
    {
        return PositionArray(0,1);
    }
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJY_() const
    {
        return PositionArray(0,1);
    }
    //we don't want these autogenerated ctors and operators
    Injection(const Injection& rhs);
    Injection& operator=(const Injection& rhs);
public:
    /**
     * \brief The constructor of a Injection object
     * 
     * Injection constructs a Injection object with:
     * \param[in] weight    the weight to do injection with (default 1.0)
     */
    Injection(const Precision weight=1.0)
        : weight_(weight),
          i_(initI_(weight)),jx_(initJX_()),jy_(initJY_()) {}

    virtual ~Injection(){}
    
    /**
     * \brief restriction() restricts the given vector to a smaler gird
     * 
     * restriction() restricts the given vector which represents a rectangular
     * grid to a smaler grid whith the double step size with weighted injection.
     * The weight is set with the constructor
     * 
     * \param[in] u                 the vector to restrict
     * \param[in] stencil           the stencil rep. of the pde needed
     *                              for matrix dep. Restrictions (not used)
     * \param[in] prolongation      Prolongation used, needed for matrix dep.
     *                              Restrictions (not used)
     * \param[in] nx                number of steps in x direction
     * \param[in] ny                number of steps in y direction
     * \throw std::domain_error     if nx or ny is not divedable by 2
     * \return                      a vector with the values on the restricted
     *                              grid
     */
    NumericArray restriction(
        const NumericArray& u,
        const Stencil& stencil,
        const Prolongation& prolongation, 
        const Index nx, const Index ny) const;
    const NumericArray& getI(
        const Position,
        const Index, 
        const Index,
        const Index,
        const Index,
        const Stencil&) const
    {
        return i_;  
    }
    const PositionArray& getJx( const Position ) const
    {
        return jx_; 
    }
    const PositionArray& getJy( const Position ) const
    {
        return jy_;
    }
};

}

#endif /*INJECTION_H_*/
