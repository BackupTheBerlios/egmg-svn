/** \file residuum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function residuum.
 * \see residuum.h
 */
#include "residuum.h"

namespace mg
{
    NumericArray residuum(
        const NumericArray& u,
        const NumericArray& fv,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny)
    {
        NumericArray result(0.0,u.size());
        if (stencil.size()<2)
            for (size_t j=1; j<ny; j++)
                for(size_t i=1; i<nx; i++)
                    result[j*(nx+1)+i]=
                            fv[j*(nx+1)+i]-stencil.apply(u,C,i,j,nx,ny);
        else
        {
            //south west corner
            result[1*(nx+1)+1]=
                    fv[1*(nx+1)+1]-stencil.apply(u,SW,1,1,nx,ny);
            //south east corner
            result[1*(nx+1)+(nx-1)]=
                    fv[1*(nx+1)+(nx-1)]-stencil.apply(u,SE,nx-1,1,nx,ny);
            //north west corner
            result[(nx-1)*(nx+1)+1]=
                    fv[(nx-1)*(nx+1)+1]-stencil.apply(u,NW,1,ny-1,nx,ny);
            //north east corner
            result[(nx-1)*(nx+1)+(nx-1)]=
                    fv[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(u,NE,nx-1,ny-1,nx,ny);
            //south boarder
            for (size_t i=2; i<nx-1; i++)
                result[1*(nx+1)+i]=
                        fv[1*(nx+1)+i]-stencil.apply(u,S,i,1,nx,ny);
            //north boarder
            for (size_t i=2; i<nx-1; i++)
                result[(nx-1)*(nx+1)+i]=
                        fv[(nx-1)*(nx+1)+i]-stencil.apply(u,N,i,ny-1,nx,ny);
            //west boarder
            for (size_t j=2; j<ny-1; j++)
                result[j*(nx+1)+1]=
                        fv[j*(nx+1)+1]-stencil.apply(u,W,1,j,nx,ny);
            //east boarder;
            for (size_t j=2; j<ny-1; j++)
                result[j*(nx+1)+(nx-1)]=
                        fv[j*(nx+1)+(nx-1)]-stencil.apply(u,E,nx-1,j,nx,ny);
            //the center
            for (size_t j=2; j<ny-1; j++)
                for (size_t i=2; i<nx-1; i++)
                    result[j*(nx+1)+i]=
                            fv[j*(nx+1)+i]-stencil.apply(u,C,i,j,nx,ny);
        }
        return result;
    }
}
