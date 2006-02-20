/** \file Stencil.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 */
#ifndef STENCIL_H_
#define STENCIL_H_

#include<valarray>
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

/**
 * \brief 	Stencil is the interface of represents a discrete differential operator.
 */
class Stencil
{
public:
	virtual ~Stencil() {}
	/**
	 * \brief applys the stencil to the given vector u at a center point
	 *
	 * apply_c evaluates \f$L_hu_h\f$ at
	 * \f$\left(i\frac{1.0}{Nx},j\frac{1}{Ny}\right)\f$
	 * The user have to take care if this point is not a point near the boarder.
	 * If it is a point near the boarder pleas use:\n
	 * E.g. apply_nw if it is the point in the north west corner.\n
	 * There is one apply for each value in the enumeration pos.
	 * \see pos
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_c(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at a west boarder point
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_w(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at the north west corner
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_nw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at a north boarder point
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_n(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at the north east corner
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_ne(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at a east boarder point
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_e(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at the south east corner
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_se(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at a south boarder point
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_s(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief applys the stencil to the given vector u at the south west corner
	 *
	 * \see apply_c
	 * 
	 * \param[in] u		the vector to apply the stencil to
	 * \param[in] i		the x position to apply the stencil to 
	 * \param[in] j		the y position to apply the stencil to
	 * \param[in] Nx	number of points in x direction (stride)
	 * \param[in] Ny	number of points in y direction
	 * \return			the value of \f$ Lu \f$
	 */
	virtual precision apply_sw(const std::valarray<precision>& u,
								const size_t i, const size_t j,
								const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at a center point
	 * 
	 * get_center returns the center element of the stencil. Because for some
	 * stencils it can make a diffrence if they are applied to a point on the
	 * boarder or the center there is one get_center for each value in the
	 * enumeration pos.
	 * \see pos
	 * \see apply_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_c(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at a west boarder point
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_w(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at the north west corner
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_nw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at a north boarder point
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_n(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at the north east corner
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_ne(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at a east boarder point
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_e(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at the south east boarder
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_se(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at a south boarder point
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_s(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coefficient of the center element at the south west corner
	 * 
	 * \see get_center_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the center coefficient
	 */
	virtual precision get_center_sw(const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at a center point
	 * 
	 * get_L_c() returns the coefficients of the stencil as valarray at a
	 * center point. L[0] is always the coefficient of the center element
	 * of the stencil the elements along the x and y axis are numbered as
	 * seen below: 
	 * \f[
	 * \left[
	 * \begin{array}{ccccccccccccccc}
	 *   & 2 &   &   & &   &   & 2 &   &   & &   &   & 2 &  \\
	 * 1 & 0 & 3 & 5 & & 5 & 1 & 0 & 3 & 6 & & 5 & 1 & 0 & 3\\
	 *   & 4 &   &   & &   &   & 4 &   &   & &   &   & 4 &  \\
	 *   & 6 &   &   & &   &   & 7 &   &   & &   &   & 6 &  \\
	 *   &   &   &   & &   &   &   &   &   & &   &   &   &  \\
	 *   & 5 &   &   & &   &   & 6 &   &   & &   &   & 6 &  \\
	 *   & 2 &   &   & &   &   & 2 &   &   & &   &   & 2 &  \\
	 * 1 & 0 & 3 & 6 & & 5 & 1 & 0 & 3 & 7 & & 5 & 1 & 0 & 3\\
	 *   & 4 &   &   & &   &   & 4 &   &   & &   &   & 4 &  \\
	 *   & 7 &   &   & &   &   & 8 &   &   & &   &   & 7 &  \\
	 *   &   &   &   & &   &   &   &   &   & &   &   &   &  \\
	 *   & 5 &   &   & &   &   & 6 &   &   & &   &   & 6 &  \\
	 *   & 2 &   &   & &   &   & 2 &   &   & &   &   & 2 &  \\
	 * 1 & 0 & 3 & 6 & & 5 & 1 & 0 & 3 & 7 & & 5 & 1 & 0 & 3\\
	 *   & 4 &   &   & &   &   & 4 &   &   & &   &   & 4 &  \\
	 * \end{array}
	 * \right]
	 * \f]
	 * (This matrix shows the numbering of all get_L return values.) 
	 * All other elements are stored in unspecific order. To get the coordinates
	 * of the other entries use get_J_x and get_J_y.\n
	 * A Stencil also have to give the garanty that the elements 0-4 are all set
	 * if size is one and that all values mentioned in the above matrix are set
	 * if the size is two.\n
	 * Because for some stencils it can make a diffrence if they are applied
	 * to a point on the boarder or the center there is one get_L for each value
	 * in the enumeration pos.\n
	 * E.g. Laplacian2D2 with the stepsize hx=hy=1 and a_x=a_y=1.0 looks like:
	 * \f[
	 * \left[\begin{array}{ccc}
	 * 		&-1	&	\\
	 * -1	& 4	&	-1\\
	 * 		&-1	&	\\
	 * \end{array}\right]
	 * \f]
	 * So we have:\n
	 * L	= {4,-1,-1,-1,-1}\n
	 * J_x 	= {0,-1,0,1,0}\n
	 * J_y 	= {0,0,1,0,-1}\n
	 * \see get_J_x
	 * \see get_J_y
	 * \see pos
	 * \see apply_c
	 * \see Laplacian2D2
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_c(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at a west boarder point
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_w(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at the north west corner
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_nw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at a north boarder point
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_n(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at the north east corner
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_ne(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at a east boarder point
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_e(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at the south east corner
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_se(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at a south boarder point
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_s(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief retruns the coefficents of the stencil at the south west corner
	 * 
	 * \see get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	the step size in x direction
	 * \param[in] Ny	the step size in y direction
	 * \return			the coefficients of Laplacian2D2
	 */
	virtual const std::valarray<precision>& get_L_sw(
							const size_t i, const size_t j,
							const size_t Nx, const size_t Ny) const =0;
	/**
	 * \brief returns the coordinate vector in x dir. for the coefficient vector
	 * \see get_L	for a more detailed description
	 * 
	 * \param[in] p		the relative position to get the coordinates vector at
	 * \return			the coordinates vector in x direction
	 */
	virtual const std::valarray<int>& get_J_x(const pos p) const =0;
	/**
	 * \brief returns the coordinate vector in y dir. for the coefficient vector
	 * \see get_L	for a more detailed description
	 * 
	 * \param[in] p		the relative position the get the coordinates vector at
	 * \return			the coordinates vector in y direction
	 */
	virtual const std::valarray<int>& get_J_y(const pos p) const =0;
	/**
	 * \brief pushs a new Prolongation on the stack of prolongations
	 * 
	 * push_prolongation is needed for galerkin like operators. Each time
	 * we are going to a coarser grid we need to tell the galerkin operator
	 * what Prolongation we are using to get back to the finer grid.
	 * 
	 * \param[in] prolong	the Prolongation to put on the stack
	 */
	virtual void push_prolongation(const Prolongation& prolong) =0;
	/**
	 * \brief removes the top Prolongation from the stack of prolongations
	 * 
	 * \see push_prolongation
	 */
	virtual void pop_prolongation() =0;
	/**
	 * \brief pushs a new restriction on the stack of restrictions
	 * 
	 * push_restriction is needed for galerkin like operators. Each time
	 * we are going to a coaser grid we need to tell the galerkin operator
	 * what Restriction we are using for this
	 * 
	 * \param[in] restriction	the Restriction to put on the stack
	 */
	virtual void push_restriction(const Restriction& restriction) =0;
	/**
	 * \brief removes the top Restriction from the stack of restrictions
	 * 
	 * \see push_restriction
	 */
	virtual void pop_restriction() =0;
	/**
	 * \brief returns the size of the stencil
	 * 
	 * The size of the stencil is the max entry of the valarrays returnd by
	 * get_J_x and get_J_y.
	 * \see get_J_x
	 * \see get_J_y
	 * 
	 * \return	the size of the stencil
	 */
	virtual size_t size() const =0;
	/**
	 * \brief is_constant says if the stencil has constant coefficients
	 * 
	 * \return		true if the stencil has constant coefficeint
	 */
	virtual bool is_constant() const =0;
};

}

#endif /*STENCIL_H_*/
