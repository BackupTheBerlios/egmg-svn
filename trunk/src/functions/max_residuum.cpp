/** \file max_residuum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function max_residuum
 * \see max_residuum.h
 */

#include <utility>
#include "max_residuum.h"

namespace mg
{
    Precision max_residuum(
        const std::valarray<Precision>& u,
        const std::valarray<Precision>& fv,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny)
    {
        Precision result=0;

        if (stencil.size()<2)
            for (size_t j=1; j<ny; j++)
                for (size_t i=1; i<nx; i++)
                {
                    Precision temp_res=
                            fv[j*(nx+1)+i]-stencil.apply(u,c,i,j,nx,ny);
                    result=std::max(std::fabs(temp_res),result);
                }
        else
        {
            //south west corner
            Precision temp_res=
                    fv[1*(nx+1)+1]-stencil.apply(u,sw,1,1,nx,ny);
            result=std::max(std::fabs(temp_res),result);
            //south east corner
            temp_res=fv[1*(nx+1)+(nx-1)]-stencil.apply(u,se,nx-1,1,nx,ny);
            result=std::max(std::fabs(temp_res),result);
            //north west corner
            temp_res=fv[(nx-1)*(nx+1)+1]-stencil.apply(u,nw,1,ny-1,nx,ny);
            result=std::max(std::fabs(temp_res),result);
            //north east corner
            temp_res=
                    fv[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(u,ne,nx-1,ny-1,nx,ny);
            result=std::max(std::fabs(temp_res),result);
            //south boarder
            for (size_t i=2; i<nx-1; i++)
            {
                temp_res=fv[1*(nx+1)+i]-stencil.apply(u,s,i,1,nx,ny);
                result=std::max(std::fabs(temp_res),result);
            }
            //north boarder
            for (size_t i=2; i<nx-1; i++)
            {
                temp_res=fv[(nx-1)*(nx+1)+i]-stencil.apply(u,n,i,ny-1,nx,ny);
                result=std::max(std::fabs(temp_res),result);
            }
            //west boarder
            for (size_t j=2; j<ny-1; j++)
            {
                temp_res=fv[j*(nx+1)+1]-stencil.apply(u,w,1,j,nx,ny);
                result=std::max(std::fabs(temp_res),result);
            }
            //east boarder
            for (size_t j=2; j<ny-1; j++)
            {
                temp_res=fv[j*(nx+1)+(nx-1)]-stencil.apply(u,e,nx-1,j,nx,ny);
                result=std::max(std::fabs(temp_res),result);
            }
            //the center
            for (size_t j=2; j<ny-1; j++)
                for (size_t i=2; i<nx-1; i++)
                {
                    temp_res=fv[j*(nx+1)+i]-stencil.apply(u,c,i,j,nx,ny);
                    result=std::max(std::fabs(temp_res),result);
                }
        }
        return result;
    }
}
