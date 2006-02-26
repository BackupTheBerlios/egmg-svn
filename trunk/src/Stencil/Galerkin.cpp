#include "Galerkin.h"

namespace mg
{
std::vector<std::valarray<int> > Galerkin::initJx_(const Stencil& sten)
{
	std::vector<std::valarray<int> > jx(9);
	for (size_t i=0; i<9; ++i)
	{
		std::valarray<int> temp = sten.getJx(static_cast<Position>(i));
		jx[i].resize(temp.size());
		jx[i] = temp;
	}
	return jx;
}

std::vector<std::valarray<int> > Galerkin::initJy_(const Stencil& sten)
{
	std::vector<std::valarray<int> > jy(9);
	for (size_t i=0; i<9; ++i)
	{
		std::valarray<int> temp = sten.getJy(static_cast<Position>(i));
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

virtual const std::valarray<Precision>& Galerkin::getL(
    const Position,
    const size_t,
    const size_t,
    const size_t,
    const size_t) const
{
    if (data_[c].find(Quadruple(i,j,Nx,Ny)) == data_[c].end())
            data_[c].insert(Quadruple(i,j,Nx,Ny),compute_L(c,i,j,Nx,Ny);
    return data_[c][Quadruple(i,j,Nx,Ny)];
}

}
