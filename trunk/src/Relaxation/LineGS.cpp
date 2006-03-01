#include "LineGS.h"

namespace mg
{
void ZebraLineGS::relax(
    std::valarray<Precision> &u,
    const std::valarray<Precision> &f, 
    const Stencil &stencil,
    const size_t nx,
    const size_t ny) const
{
    // valarrays needed for LR-decomposition of a tridiagonal matrix
    std::valarray<Precision> rhs(0.0,u.size());
    switch (stencil.size())
    {
        case 1:  // stencil of size 1
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    ninepointxline(u, f, rhs, stencil, nx, ny);
                    ninepointyline(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case XDIR:
                {
                    ninepointxline(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case YDIR:
                {
                    ninepointyline(u, f, rhs, stencil, nx, ny);
                    break;
                }
                default:
                {
                    std::cerr << "Error in direction of the line relaxation!\n";
                    break;
                }
            }
            break;          
        }
        case 2:  //stencil of size 2
        {
            switch (direction_)
            {
                case ALTDIR:
                { 
                    xline(u, f, rhs, stencil, nx, ny);
                    yline(u, f, rhs, stencil, nx, ny);
                    break;
                }
                case XDIR:
                {
                    xline(u, f, rhs, stencil, nx, ny);
                    break;
                }   
                case YDIR:
                {
                    yline(u, f, rhs, stencil, nx, ny);
                    break;
                }
                default:
                {
                    std::cerr << "Error in direction of the line relaxation!\n";
                    break;
                }
            }
            break;
        }
        default:
        {
            std::cerr << "Stencil is too big (size>2)!" << std::endl;
            break;
        }
    }
}
}
