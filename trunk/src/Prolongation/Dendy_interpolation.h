/** \file Dendy_interpolation.h
 * \author Benedikt Engbroks
 * \brief Contains the interface of the class Dendy_interpolation
 * 
 * This file contains the interface of Dendy_interpolation. The
 * implementation is in Dendy_interpolation.cpp
 */
#ifndef DENDY_INTERPOLATION_H_
#define DENDY_INTERPOLATION_H_

#include "Prolongation.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief Dendy_interpolation is a matrix-dependent 2D Prolongation Operator
 * 
 * Dendy_interpolation represents a 2D Prolongation Operator that uses
 * matrix-dependent interpolation by Dendy to do its job.
 */
class Dendy_interpolation : public mg::Prolongation
{
private:
	std::valarray<precision> t;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	//initialize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		const int t[] = {0,-1,0,1,0,-1,1,1,-1};
		return std::valarray<int>(t,9);
	}
	//initialize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		const int t[] = {0,0,1,0,-1,1,1,-1,-1};
		return std::valarray<int>(t,9);
	}
	//we don't want these autogenerated contrs and operators
	Dendy_interpolation(const Dendy_interpolation& rhs);
	Dendy_interpolation& operator=(const Dendy_interpolation& rhs);

	static const int SW = 0;
	static const int S = 1;
	static const int SE = 2;
	static const int W = 3;
	static const int C = 4;
	static const int E = 5;
	static const int NW = 6;
	static const int N = 7;
	static const int NE = 8;

public:                                                              
	Dendy_interpolation() : t(9), J_x(init_J_x()), J_y(init_J_y()) {}
	virtual ~Dendy_interpolation() {}
	/**
	 * \brief prolong does a matrix-dependent interpolation on the input vector
	 * 
	 * \param u		the vector representing a rectangle to prolongate
	 * \param Nx	Number of steps in x direction
	 * \param Ny	Number of steps in y direction
	 * \return 		a vector representing the prolongated rectangle of 
	 * 				size 2*(Nx+1)*2*(Ny+1)
	 */
	std::valarray<precision> prolong(const std::valarray<precision>& u, const Stencil& stencil, 
									const size_t Nx,const size_t Ny) const;
									
	const std::valarray<precision>& get_I(const size_t i, const size_t j, 
		const size_t Nx, const size_t Ny, const Stencil& stencil)
	{
		std::valarray<precision> L = stencil.get_L_c(0,0,Nx,Ny);
		std::valarray<int> J_x = stencil.get_J_x(c);
		std::valarray<int> J_y = stencil.get_J_y(c);
		std::valarray<size_t> posi(9);
		for (size_t jj=0;jj<J_x.size();jj++)
		{
			posi[(J_x[jj]+1)+3*(J_y[jj]+1)] = jj;
		}

		precision scale = 0;
		precision weight1 = 0;
		precision weight2 = 0;
		precision erg = 0;

		// C
		t[0] = 1.0;
		// W
			L = stencil.get_L_c(i-1,j,Nx,Ny);
			scale = -L[0];
			weight2 = 0;
			if (posi[S] != 0) scale -= L[posi[S]];
			if (posi[N] != 0) scale -= L[posi[N]];
			if (posi[E] != 0) weight2 += L[posi[E]];
			if (posi[SE] != 0) weight2 += L[posi[SE]];
			if (posi[NE] != 0) weight2 += L[posi[NE]];
			t[1] = weight2 / scale; 
		// N
			L = stencil.get_L_c(i,j+1,Nx,Ny);
			scale = -L[0];
			weight1 = 0;
			if (posi[W] != 0) scale -= L[posi[W]];
			if (posi[E] != 0) scale -= L[posi[E]];
			if (posi[S] != 0) weight1 += L[posi[S]];
			if (posi[SW] != 0) weight1 += L[posi[SW]];
			if (posi[SE] != 0) weight1 += L[posi[SE]];
			t[2] = weight1 / scale;
		// E
			L = stencil.get_L_c(i+1,j,Nx,Ny);
			scale = -L[0];
			weight1 = 0;
			if (posi[S] != 0) scale -= L[posi[S]];
			if (posi[N] != 0) scale -= L[posi[N]];
			if (posi[W] != 0) weight1 += L[posi[W]];
			if (posi[SW] != 0) weight1 += L[posi[SW]];
			if (posi[NW] != 0) weight1 += L[posi[NW]];
			t[3] = weight1 / scale; 		
		// S
			L = stencil.get_L_c(i,j-1,Nx,Ny);
			scale = -L[0];
			weight2 = 0;
			if (posi[W] != 0) scale -= L[posi[W]];
			if (posi[E] != 0) scale -= L[posi[E]];
			if (posi[N] != 0) weight2 += L[posi[N]];
			if (posi[NW] != 0) weight2 += L[posi[NW]];
			if (posi[NE] != 0) weight2 += L[posi[NE]];
			t[4] = weight2 / scale;
		// NW
			L = stencil.get_L_c(i-1,j+1,Nx,Ny);
			scale = -L[0];
			erg = 0;
			if (posi[E] != 0) erg += L[posi[E]] * t[2];
			if (posi[S] != 0) erg += L[posi[S]] * t[1];
			if (posi[SE] != 0) erg += L[posi[SE]];
			t[5] = erg / scale;
		// NE
			L = stencil.get_L_c(i+1,j+1,Nx,Ny);
			scale = -L[0];
			erg = 0;
			if (posi[W] != 0) erg += L[posi[W]] * t[2];
			if (posi[S] != 0) erg += L[posi[S]] * t[3];
			if (posi[SW] != 0) erg += L[posi[SW]];
			t[6] = erg / scale;
		// SE
			L = stencil.get_L_c(i+1,j-1,Nx,Ny);
			scale = -L[0];
			erg = 0;
			if (posi[W] != 0) erg += L[posi[W]] * t[4];
			if (posi[N] != 0) erg += L[posi[N]] * t[3];
			if (posi[NW] != 0) erg += L[posi[NW]];
			t[7] = erg / scale;
		// SW
			L = stencil.get_L_c(i-1,j-1,Nx,Ny);
			scale = -L[0];
			erg = 0;
			if (posi[E] != 0) erg += L[posi[E]] * t[4];
			if (posi[N] != 0) erg += L[posi[N]] * t[1];
			if (posi[NE] != 0) erg += L[posi[NE]];
			t[8] = erg / scale;		
		return t;			
	}
	const std::valarray<int>& get_J_x() const
	{
		return J_x;	
	}
	const std::valarray<int>& get_J_y() const
	{
		return J_y;
	}
};

}

#endif /*DENDY_INTERPOLATION_H_*/
