/** \file Transfer_Operator.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the abstract class Transfer_Operator
 */
#ifndef TRANSFER_OPERATOR_H_
#define TRANSFER_OPERATOR_H_

namespace mg
{
class Stencil;
/**
 * \brief Transfer_Operator is the interface of a 2D transfer operator
 */
class Transfer_Operator
{
	/**
	 * \brief get_I returns the stencil representation of the Transfer_Operator
	 * 
	 * Retruns the coefficents of the stencil representation as valarray. For a
	 * 9 Point Stencil the entries are numbered like this:
	 * \f[
	 * \left[
	 * \begin{array}{ccc}
	 * 5 & 2 & 6 \\
	 * 1 & 0 & 3 \\
	 * 8 & 4 & 7 \\
	 * \end{array}
	 * \right]
	 * =
	 * \left[
	 * \begin{array}{ccc}
	 * nw & n & ne \\
	 * w & c & e \\
	 * sw & s & se \\
	 * \end{array}
	 * \right]
	 * \f]
	 * The enum pos has the same ordering, please use the enum values for
	 * addresing. It is guaranteed that all 9 Values are set. To iterate over
	 * the entries you can also use the valarrays returnd by get_J_x and
	 * get_J_y. This is analog to the get_L and get_J_x get_J_y Methods of the
	 * Stencil interface.
	 * 
	 * \see pos
	 * \see Stencil.get_L_c
	 * 
	 * \param[in] i		the x position
	 * \param[in] j		the y position
	 * \param[in] Nx	number of points in x direction
	 * \param[in] Ny	number of points in y direction
	 * \return			the coefficients of the Transfer_Operator
	 */
	virtual const std::valarray<precision>& get_I(
									const size_t i, const size_t j,
									const size_t Nx, const size_t Ny, 
									const Stencil& stencil) const =0;
	/**
	 * \brief returns the coordinate vector in x dir. for the coefficient vec. I
	 * \see get_I
	 * 
	 * \return			the coordinates vector in y direction
	 */
	virtual const std::valarray<int>& get_J_x() const =0;
	/**
	 * \brief returns the coordinate vector in y dir. for the coefficient vec. I
	 * \see get_I
	 * 
	 * \return			the coordinates vector in y direction
	 */
	virtual const std::valarray<int>& get_J_y() const =0;
};

}

#endif /*TRANSFER_OPERATOR_H_*/
