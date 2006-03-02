#include "Galerkin.h"

namespace mg
{
std::vector<PositionArray > Galerkin::initJx_(const Stencil& sten)
{
	std::vector<PositionArray > jx(9);
	for (size_t i=0; i<9; ++i)
	{
		PositionArray temp = sten.getJx(static_cast<Position>(i));
		jx[i].resize(temp.size());
		jx[i] = temp;
	}
	return jx;
}

std::vector<PositionArray > Galerkin::initJy_(const Stencil& sten)
{
	std::vector<PositionArray > jy(9);
	for (size_t i=0; i<9; ++i)
	{
		PositionArray temp = sten.getJy(static_cast<Position>(i));
		jy[i].resize(temp.size());
		jy[i] = temp;
	}
	return jy;
}

void Galerkin::updateSize_()
{
	
}
void Galerkin::updateJxJy_()
{
	
}

virtual const NumericArray& Galerkin::getL(
    const Position,
    const size_t,
    const size_t,
    const size_t,
    const size_t) const
{
    if (data_[C].find(Quadruple(i,j,Nx,Ny)) == data_[C].end())
            data_[C].insert(Quadruple(i,j,Nx,Ny),compute_L(C,i,j,Nx,Ny);
    return data_[C][Quadruple(i,j,Nx,Ny)];
}

}
