#include "Galerkin.h"

namespace mg
{
std::vector<std::valarray<int> > Galerkin::init_J_x(const Stencil& sten)
{
	std::vector<std::valarray<int> > Jx(9);
	for (size_t i=0;i<9;i++)
	{
		std::valarray<int> temp = sten.getJx(static_cast<Position>(i));
		Jx[i].resize(temp.size());
		Jx[i] = temp;
	}
	return Jx;
}

std::vector<std::valarray<int> > Galerkin::init_J_y(const Stencil& sten)
{
	std::vector<std::valarray<int> > Jy(9);
	for (size_t i=0;i<9;i++)
	{
		std::valarray<int> temp = sten.getJy(static_cast<Position>(i));
		Jy[i].resize(temp.size());
		Jy[i] = temp;
	}
	return Jy;
}

void Galerkin::update_size()
{
	
}
void Galerkin::update_J_x_J_y()
{
	
}
}
