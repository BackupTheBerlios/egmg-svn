/** \file directSolver.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief directSolver.cpp contains the impl. of the func. directSolver.
 */

#include <algorithm>
#include "directSolver.h"
#include "../Stencil/Stencil.h"
#include "../general/parameters.h"

namespace mg
{

namespace
{

PositionArray pivotLU(
    NumericArray& matrix,
    const Index dimension)
{
    PositionArray permutation(dimension);
    for (Index sx = 0; sx < dimension; ++sx)
        permutation[sx]=static_cast<Integer>( sx );

    for (Index sx=0; sx<dimension; ++sx)
    {
        Index max_zeile=sx;

        for (Index sy=sx; sy<dimension; ++sy)
            if (std::abs(matrix[permutation[sy       ]*dimension+sx])
               >std::abs(matrix[permutation[max_zeile]*dimension+sx]))
                max_zeile=sy;

        std::swap(permutation[sx],permutation[max_zeile]);

        for (Index sy=sx+1; sy<dimension; ++sy)
        {
            matrix[permutation[sy]*dimension+sx]/=
                    matrix[permutation[sx]*dimension+sx];
            Precision pivot=matrix[permutation[sy]*dimension+sx];

            for (Index k=sx+1; k<dimension; k++)
                matrix[permutation[sy]*dimension+k]-=
                     pivot
                    *matrix[permutation[sx]*dimension+k];
        }
    }
    return permutation;
}

NumericArray solve(
    NumericArray matrix,
    const NumericArray& b)
{
    const Index dimension = b.size();
    PositionArray permutation=pivotLU(matrix,dimension);

    NumericArray result(dimension);
    //solve Ly=b
    for (Index sx=0; sx<dimension; ++sx)
    {
        result[sx]=b[permutation[sx]];
        for (Index sy=0; sy<sx; ++sy)
            result[sx]-=matrix[permutation[sx]*dimension+sy]*result[sy];
    }
    //solve Rx=y
    for (Index sx=0; sx<dimension; ++sx)
    {
        for (Index sy=0; sy<sx; ++sy)
            result[dimension-1-sx]-=
                 matrix[permutation[dimension-1-sx]*dimension+dimension-1-sy]
                *result[dimension-1-sy];
        result[dimension-1-sx]/=
            matrix[(permutation[dimension-1-sx]+1)*dimension-1-sx];
    }
    return result;
}

}

NumericArray directSolver(
    const LinearEquationSystem& les)
{  
    return solve(les.first,les.second);
}
  
}
