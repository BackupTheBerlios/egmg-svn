/** \file Helmholtz2D2.h
 * \author Andr� Oeckerath
 * \brief Contains the class Helmholtz2D2
 */
#ifndef Helmholtz2D2_H_
#define Helmholtz2D2_H_

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

/**
 * \brief 	Helmholtz2D2 represents a discrete laplace like operator of second
 * 			order.
 * 
 * Helmholtz2D2 is the stencil representing the discrete differential operator
 * \f[
 *	L_h u_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}+c_{xy}(u_h)
 * \f]
 */
class Helmholtz2D2 : public mg::Stencil
{
private:
	mutable std::valarray<Precision> L;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	const Precision ax;
	const Precision ay;
	const Precision c_;
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
	//we don't want these autogenerated ctors and operators
	Helmholtz2D2(const Helmholtz2D2& rhs);
	Helmholtz2D2& operator=(const Helmholtz2D2& rhs);
public:
	/**
	 * \brief The constructor of a Helmholtz2D2 object
	 * 
	 * Helmholtz2D2 constructs a Helmholtz2D2 where \f$a_x\f$ , \f$a_y\f$
	 * and \f$c_{xy}\f$ are given by:
	 * \param[in] a_x	coefficient of the diff. operator (default 1.0)
	 * \param[in] a_y	coefficient of the diff. operator (default 1.0)
	 * \param[in] c coefficient of the diff. operator (default 1.0)
	 */
	explicit Helmholtz2D2(Precision a_x =1.0,Precision a_y =1.0, Precision c_x_y = 1.0) 
		: L(5), J_x(init_J_x()), J_y(init_J_y()), ax(a_x), ay(a_y), c_(c_x_y) {}
	virtual ~Helmholtz2D2() {}
	inline Precision apply_c(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return (2.0*ax*Nx*Nx+2.0*ay*Ny*Ny+c_)*u[j*(Nx+1)+i]
				-1.0*ax*Nx*Nx*u[j*(Nx+1)+i-1]-1.0*ax*Nx*Nx*u[j*(Nx+1)+i+1]
				-1.0*ay*Ny*Ny*u[(j-1)*(Nx+1)+i]-1.0*ay*Ny*Ny*u[(j+1)*(Nx+1)+i];
	}
	inline Precision get_center_c(const size_t, const size_t,
							const size_t Nx, const size_t Ny) const
	{
		return (2.0*ax*Nx*Nx+2.0*ay*Ny*Ny+c_);
	}
	inline Precision apply_w(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_nw(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_n(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_ne(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_e(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_se(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_s(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision apply_sw(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return apply_c(u,i,j,Nx,Ny);
	}
	inline Precision get_center_w(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_nw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_n(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_ne(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_e(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_se(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_s(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline Precision get_center_sw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_center_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_c(
							const size_t, const size_t,
							const size_t Nx, const size_t Ny) const
	{
		L[0] = 2.0*ax*Nx*Nx+2.0*ay*Ny*Ny+c_;
		L[1] = L[3] = -1.0*ax*Nx*Nx;
		L[2] = L[4] = -1.0*ay*Ny*Ny;
		return L;
	}
	inline const std::valarray<Precision>& get_L_w(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_nw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_n(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_ne(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_e(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_se(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_s(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<Precision>& get_L_sw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const
	{
		return get_L_c(i,j,Nx,Ny);
	}
	inline const std::valarray<int>& get_J_x(const Position =c) const
	{
		return J_x;
	}
	inline const std::valarray<int>& get_J_y(const Position =c) const
	{
		return J_y;
	}
	/**
	 * \brief does nothing for Helmholtz2D2
	 * \see Stencil
	 */
	void push_prolongation(const Prolongation&) {}
	/**
	 * \brief does nothing for Helmholtz2D2
	 * \see Stencil
	 */
	void pop_prolongation() {}
	/**
	 * \brief does nothing for Helmholtz2D2
	 * \see Stencil
	 */
	void push_restriction(const Restriction&) {}
	/**
	 * \brief does nothing for Helmholtz2D2
	 * \see Stencil
	 */
	void pop_restriction() {}
	/**
	 * \brief gives the max expansion of Helmholtz2D2
	 * 
	 * \return	1
	 */
	inline size_t size() const
	{
		return 1;
	}
	/**
	 * \brief returns true, because Helmholtz2D2 is constant
	 * 
	 * \return	true
	 */
	inline bool is_constant() const
	{
		return true;
	}
	
};

}

#endif /*Helmholtz2D2_H_*/
