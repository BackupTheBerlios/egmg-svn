/** \file Bicubic_interpolation.h
 * \brief Contains the interface of the class Bicubic_interpolation
 * \author Benedikt Engbroks
 * 
 * This file contains the interface of Bicubic_interpolation. The
 * implementation is in Bicubic_interpolation.cpp
 */
#ifndef BICUBIC_INTERPOLATION_H_
#define BICUBIC_INTERPOLATION_H_

#include "Prolongation.h"

namespace mg
{

/**
 * \brief Bicubic_interpolation is a 2D Prolongation Operator
 * 
 * Bicubic_interpolation represents a 2D Prolongation Operator that uses
 * bicubic interpolation to do its job.
 */
class Bicubic_interpolation : public mg::Prolongation
{
private:
	const std::valarray<precision> I;
	std::valarray<precision> IRand;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	//initialize I, makes it possible to make I const
	std::valarray<precision> init_I() const
	{	// use the saved I only for inner points
		const precision t[] = {1.0,9.0/16,9.0/16,9.0/16,9.0/16,81.0/256,
			81.0/256,81.0/256,81.0/256,-1.0/16,-9.0/256,1.0/256,
			-9.0/256,-1.0/16,-9.0/256,1.0/256,-9.0/256,-1.0/16,
			-9.0/256,1.0/256,-9.0/256,-1.0/16,-9.0/256,1.0/256,
			-9.0/256};
		return std::valarray<precision>(t,25);
	}
	//initialize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		const int t[] = {0,-1,0,1,0,-1,1,1,-1,-3,-3,-3,-1,0,1,3,3,3,
				 3,3,1,0,-1,-3,-3};
		return std::valarray<int>(t,25);
	}
	//initialize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		const int t[] = {0,0,1,0,-1,1,1,-1,-1,0,1,3,3,3,3,3,1,0,-1,-3,
				 -3,-3,-3,-3,-1};
		return std::valarray<int>(t,25);
	}
	//we don't want these autogenerated contors and operators
	Bicubic_interpolation(const Bicubic_interpolation& rhs);
	Bicubic_interpolation& operator=(const Bicubic_interpolation& rhs);
public:
	Bicubic_interpolation() : I(init_I()), IRand(25), J_x(init_J_x()), J_y(init_J_y()) {}
	virtual ~Bicubic_interpolation() {}
	/**
	 * \brief prolong does a bicubic interpolation on the input vector
	 * 
	 * \param u		the vector representing a rectangel to prolongate
	 * \param Nx	Number of steps in x direction
	 * \param Ny	Number of steps in y direction
	 * \return 		a vector representing the prolongated rectangel of 
	 * 				size 2*(Nx+1)*2*(Ny+1)
	 */
	std::valarray<precision> prolong(const std::valarray<precision>& u,
									const Stencil& stencil,
									const size_t Nx,const size_t Ny) const;
	const std::valarray<precision>& get_I(const size_t i, const size_t j,
					      const size_t Nx, const size_t Ny, const Stencil&)
	{
		if (j == 1)        // unterer Rand
		{
			if (i == 1)   // untere linke Ecke
			{
				const precision t[] = {1.0,3.0/8,3.0/4,3.0/4,3.0/8,9.0/32,
					9.0/16,9.0/32,9.0/64,0.0,0.0,0.0,
					-3.0/64,-1.0/8,-3.0/32,1.0/64,-3.0/32,-1.0/8,
					-3.0/64,0.0,0.0,0.0,0.0,0.0,0.0};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
			else if (i == Nx-1)   // untere rechte Ecke
			{
				const precision t[] = {1.0,3.0/4,3.0/4,3.0/8,3.0/8,9.0/16,
					9.0/32,9.0/64,9.0/32,-1.0/8,-3.0/32,1.0/64,
					-3.0/32,-1.0/8,-3.0/64,0.0,0.0,0.0,
					0.0,0.0,0.0,0.0,0.0,0.0,-3.0/64};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
			else
			{
				const precision t[] = {1.0,9.0/16,3.0/4,9.0/16,3.0/8,27.0/64,
					27.0/64,27.0/128,27.0/128,-1.0/16,-3.0/64,1.0/128,
					-9.0/128,-1.0/8,-9.0/128,1.0/128,-3.0/64,-1.0/16,
					-3.0/128,0.0,0.0,0.0,0.0,0.0,-3.0/128};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
		}
		else if (j == Ny-1)  // oberer Rand
		{
			if (i == 1)   // obere linke Ecke
			{
				const precision t[] = {1.0,3.0/8,3.0/8,3.0/4,3.0/4,9.0/64,
					9.0/32,9.0/16,9.0/32,0.0,0.0,0.0,
					0.0,0.0,0.0,0.0,-3.0/64,-1.0/8,
					-3.0/32,1.0/64,-3.0/32,-1.0/8,-3.0/64,0.0,0.0};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
			else if (i == Nx-1)   // obere rechte Ecke
			{
				const precision t[] = {1.0,3.0/4,3.0/8,3.0/8,3.0/4,9.0/32,
					9.0/64,9.0/32,9.0/16,-1.0/8,-3.0/64,0.0,
					0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-3.0/64,-1.0/8,
					-3.0/32,1.0/64,-3.0/32};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
			else
			{
				const precision t[] = {1.0,9.0/16,3.0/8,9.0/16,3.0/4,27.0/128,
					27.0/128,27.0/64,27.0/64,-1.0/16,-3.0/128,0.0,
					0.0,0.0,0.0,0.0,-3.0/128,-1.0/16,-3.0/64,1.0/128,
					-9.0/128,-1.0/8,-9.0/128,1.0/128,-3.0/64};
				IRand = std::valarray<precision>(t,25);
				return IRand;				
			}
		}
		else if (i == 1)   // linker Rand
		{
			const precision t[] = {1.0,3.0/8,9.0/16,3.0/4,9.0/16,27.0/128,
				27.0/64,27.0/64,27.0/128,0.0,0.0,0.0,-3.0/128,-1.0/16,
				-3.0/64,1.0/128,-9.0/128,-1.0/8,-9.0/128,1.0/128,
				-3.0/64,-1.0/16,-3.0/128,0.0,0.0};
			IRand = std::valarray<precision>(t,25);
			return IRand;				
		}
		else if (i == Nx-1)  // rechter Rand
		{
			const precision t[] = {1.0,3.0/4,9.0/16,3.0/8,9.0/16,27.0/64,
				27.0/128,27.0/128,27.0/64,-1.0/8,-9.0/128,1.0/128,
				-3.0/64,-1.0/16,-3.0/128,0.0,0.0,0.0,0.0,0.0,
				-3.0/128,-1.0/16,-3.0/64,1.0/128,-9.0/128};
			IRand = std::valarray<precision>(t,25);
			return IRand;				
		}
		else
			return I;
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

#endif /*BICUBIC_INTERPOLATION_H_*/
