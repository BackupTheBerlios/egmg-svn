/** \file ConDiff2D1.h
 * \author Andr� Oeckerath
 * \brief Contains the class ConDiff2D1
 */
#ifndef ConDiff2D1_H_
#define ConDiff2D1_H_

#include <cmath>
#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"


namespace mg
{

/**
 * \brief 	ConDiff2D1 represents a discrete Convection Diffusion operator of first
 * 			order.
 * 
 * ConDiff2D1 is the stencil representing the discrete differential operator
 * \f[
 *	L_hu_h := \epsilon (-(u_h)_{xx}-(u_h)_{yy}) + cos(\beta) u_x + sin(\beta) u_y
 * \f]
 */
class ConDiff2D1 : public mg::Stencil
{
private:
	mutable std::valarray<precision> L;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	const precision epsilon;
	const precision beta;
	precision a1;
	precision a2;
	//initilize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		const int t[] = {0,-1,0,1,0};
		return std::valarray<int>(t,5);
	}
	//initilize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		const int t[] = {0,0,1,0,-1};
		return std::valarray<int>(t,5);
	}
    //initilize a1 and a2
	precision init_coeff(precision beta_) 
	{
		a1 = cos(beta_);
		a2 = sin(beta_);
		return beta_;
	}
	//we don't want these autogenerated ctors and operators
	ConDiff2D1(const ConDiff2D1& rhs);
	ConDiff2D1& operator=(const ConDiff2D1& rhs);
public:
	/**
	 * \brief The constructor of a ConDiff2D1 object
	 * 
	 * ConDiff2D1 constructs a ConDiff2D1 where \f$\epsilon\f$, \f$a_1\f$ and \f$a_2\f$
	 * are given by:
     * \param[in] \epsilon coefficient of the diffusion part (default 1.0)
	 * \param[in] a_1 = cos(\beta)	coefficient of the convection part (default \pi/4)
	 * \param[in] a_2 = sin(\beta)	coefficient of the convection part (default \pi/4)
	 */
	explicit ConDiff2D1(precision epsilon_ =1.0, precision beta_ = 0.78539816339745) 
		: L(5), J_x(init_J_x()), J_y(init_J_y()), epsilon(epsilon_), beta(init_coeff(beta_)) {}
	virtual ~ConDiff2D1() {}
	inline precision apply_c(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return (2.0*epsilon*Nx*Nx+2.0*epsilon*Ny*Ny+fabs(a1)*Nx+fabs(a2)*Ny)*u[j*(Nx+1)+i]
				+(-1.0*epsilon*Nx*Nx+(-a1-fabs(a1))*Nx/2)*u[j*(Nx+1)+i-1]
				+(-1.0*epsilon*Nx*Nx+(a1-fabs(a1))*Nx/2)*u[j*(Nx+1)+i+1]
				+(-1.0*epsilon*Ny*Ny+(-a2-fabs(a2))*Ny/2)*u[(j-1)*(Nx+1)+i]
				+(-1.0*epsilon*Ny*Ny+(a2-fabs(a2))*Ny/2)*u[(j+1)*(Nx+1)+i];
	}
	inline precision get_center_c(const size_t, const size_t,
							const size_t Nx, const size_t Ny) const
	{
		return (2.0*epsilon*Nx*Nx+2.0*epsilon*Ny*Ny+fabs(a1)*Nx+fabs(a2)*Ny);
	}
	inline precision apply_w(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_nw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_n(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_ne(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_e(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_se(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_s(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision apply_sw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline precision get_center_w(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_nw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_n(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_ne(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_e(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_se(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_s(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline precision get_center_sw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_c(
							const size_t, const size_t,
							const size_t Nx, const size_t Ny) const
	{
		L[0] = 2.0*epsilon*Nx*Nx+2.0*epsilon*Ny*Ny+fabs(a1)*Nx+fabs(a2)*Ny;
		L[1] = -1.0*epsilon*Nx*Nx+(-a1-fabs(a1))*Nx/2;
		L[3] = -1.0*epsilon*Nx*Nx+(a1-fabs(a1))*Nx/2;
		L[2] = -1.0*epsilon*Ny*Ny+(-a2-fabs(a2))*Ny/2;
		L[4] = -1.0*epsilon*Ny*Ny+(a2-fabs(a2))*Ny/2;
		return L;
	}
	inline const std::valarray<precision>& get_L_w(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_nw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_n(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_ne(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_e(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_se(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_s(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<precision>& get_L_sw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<int>& get_J_x(const pos =c) const
	{
		return J_x;
	}
	inline const std::valarray<int>& get_J_y(const pos =c) const
	{
		return J_y;
	}
	/**
	 * \brief does nothing for ConDiff2D1
	 * \see Stencil
	 */
	void push_prolongation(const Prolongation&) {}
	/**
	 * \brief does nothing for ConDiff2D1
	 * \see Stencil
	 */
	void pop_prolongation() {}
	/**
	 * \brief does nothing for ConDiff2D1
	 * \see Stencil
	 */
	void push_restriction(const Restriction&) {}
	/**
	 * \brief does nothing for ConDiff2D1
	 * \see Stencil
	 */
	void pop_restriction() {}
	/**
	 * \brief gives the max expansion of ConDiff2D1
	 * 
	 * \return	1
	 */
	inline size_t size() const
	{
		return 1;
	}
	/**
	 * \brief returns true, because ConDiff2D1 is constant
	 * 
	 * \return	true
	 */
	inline bool is_constant() const
	{
		return true;
	}
	
};

}

#endif /*ConDiff2D1_H_*/
