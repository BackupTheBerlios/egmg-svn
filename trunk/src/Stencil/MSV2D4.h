/** \file MSV2D4.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the class MSV2D4
 */
#ifndef MSV2D4_H_
#define MSV2D4_H_

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

/**
 * \brief 	MSV2D4 represents a discrete laplace like operator of second
 * 			order.
 * 
 * MSV2D4 is the stencil representing the discrete differential operator
 * \f[
 *	L_h u_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}
 * \f]
 */
class MSV2D4 : public mg::Stencil
{
private:
	mutable std::valarray<Precision> L;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	const Precision ax;
	const Precision ay;
	//initilize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		const int t[] = {0,-1,0,1,0,-1,1,-1,1};
		return std::valarray<int>(t,9);
	}
	//initilize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		const int t[] = {0,0,1,0,-1,-1,-1,1,1};
		return std::valarray<int>(t,9);
	}
	//we don't want these autogenerated ctors and operators
	MSV2D4(const MSV2D4& rhs);
	MSV2D4& operator=(const MSV2D4& rhs);
public:
	/**
	 * \brief The constructor of a MSV2D4 object
	 * 
	 * MSV2D4 constructs a MSV2D4 where \f$a_x\f$ and \f$a_y\f$
	 * are given by:
	 * \param[in] a_x	coefficient of the diff. operator (default 1.0)
	 * \param[in] a_y	coefficient of the diff. operator (default 1.0)
	 */
	explicit MSV2D4(Precision a_x =1.0,Precision a_y =1.0) : L(9), J_x(init_J_x()), J_y(init_J_y()), ax(a_x), ay(a_y) {}
	virtual ~MSV2D4() {}
	inline Precision apply_c(const std::valarray<Precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const
	{
		return ((20.0*Nx*Nx+20.0*Ny*Ny)/12.0)*u[j*(Nx+1)+i]
				-10.0*Nx*Nx/12.0*u[j*(Nx+1)+i-1]-10.0*Nx*Nx/12.0*u[j*(Nx+1)+i+1]
				-10.0*Ny*Ny/12.0*u[j*(Nx+1)+i-(Nx+1)]-10.0*Ny*Ny/12.0*u[j*(Nx+1)+i+(Nx+1)]
				+2.0*Nx*Nx/12.0*u[j*(Nx+1)+i-(Nx+1)]+2.0*Nx*Nx/12.0*u[j*(Nx+1)+i+(Nx+1)]
				+2.0*Ny*Ny/12.0*u[j*(Nx+1)+i-1]+2.0*Ny*Ny/12.0*u[j*(Nx+1)+i+1]
				-1.0*Nx*Nx/12.0*(u[j*(Nx+1)+i-(Nx+1)-1] + u[j*(Nx+1)+i-(Nx+1)+1]
				                  +u[j*(Nx+1)+i+(Nx+1)-1] + u[j*(Nx+1)+i+(Nx+1)+1])
				-1.0*Ny*Ny/12.0*(u[j*(Nx+1)+i-(Nx+1)-1] + u[j*(Nx+1)+i-(Nx+1)+1]
				                  +u[j*(Nx+1)+i+(Nx+1)-1] + u[j*(Nx+1)+i+(Nx+1)+1]);
	}
	inline Precision get_center_c(const size_t, const size_t,
							const size_t Nx, const size_t Ny) const
	{
		return ((20.0*Nx*Nx+20.0*Ny*Ny)/12.0);
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
		L[0] = 20.0*Nx*Nx/12+20.0*Ny*Ny/12;
		L[1] = L[3] = -10.0*Nx*Nx/12 + 2.0*Ny*Ny/12;
		L[2] = L[4] = -10.0*Ny*Ny/12 + 2.0*Nx*Nx/12;
		L[5] = L[6] = L[7] = L[8] = -1.0*Nx*Nx/12 - 1.0*Ny*Ny/12;
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
	 * \brief does nothing for MSV2D4
	 * \see Stencil
	 */
	void push_prolongation(const Prolongation&) {}
	/**
	 * \brief does nothing for MSV2D4
	 * \see Stencil
	 */
	void pop_prolongation() {}
	/**
	 * \brief does nothing for MSV2D4
	 * \see Stencil
	 */
	void push_restriction(const Restriction&) {}
	/**
	 * \brief does nothing for MSV2D4
	 * \see Stencil
	 */
	void pop_restriction() {}
	/**
	 * \brief gives the max expansion of MSV2D4
	 * 
	 * \return	1
	 */
	inline size_t size() const
	{
		return 1;
	}
	/**
	 * \brief returns true, because MSV2D4 is constant
	 * 
	 * \return	true
	 */
	inline bool is_constant() const
	{
		return true;
	}
	
};

}

#endif /*LAPLACIAN2D2_H_*/
