/** \file Galerkin.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface the class Galerkin.
 */
#ifndef GALERKIN_H_
#define GALERKIN_H_

#include "PreCalculatedStencil.h"
#include <vector>
#include <ostream>

namespace mg
{

class Galerkin : public PreCalculatedStencil
{  
public:
    Galerkin( const Stencil& fineGridOperator );
    
    virtual ~Galerkin();

    
protected:
    virtual void update(
        const Restriction& restriction,
        const Prolongation& prolongation,
        const Index nx,
        const Index ny );
    virtual void update();

private:
    std::vector<const Prolongation*> prolongations_;
    std::vector<const Restriction*> restrictions_;
};

bool CheckSymmetry(const Restriction& restriction, Galerkin& stencil, const Prolongation& prolongation, const Index nx, const Index ny, std::ostream& out );

}

#endif /*GALERKIN_H_*/
