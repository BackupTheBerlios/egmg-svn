/** \file Stencil.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Stencil.
 */
#ifndef STENCIL_H_
#define STENCIL_H_

#include "../general/parameters.h"
#include "../general/DiscreteFunction.h"

namespace mg
{
    
class Prolongation;
class Restriction;
class DiscreteFunction;

/**
 * \brief   Stencil is the interface of represents a discrete differential operator.
 */
class Stencil
{
public:
	Stencil();
    virtual ~Stencil();
    
    /**
     * \brief applys the stencil to the given vector u at the given point
     *
     * apply evaluates \f$L_h u_h\f$ at
     * \f$\left(sx\frac{1.0}{nx},sy\frac{1}{ny}\right)\f$
     * The user have to take care if this point is a point near the boarder and
     * set the Position parameter pos correctly.
     * If it is a point near the boarder pleas use:\n
     * E.g. apply(u,nw,sx,sy,nx,ny) if it is the point in the north west corner.\n
     * \see Position
     * 
     * \param[in] u         the vector to apply the stencil to
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position to apply the stencil to 
     * \param[in] sy        the y position to apply the stencil to
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the value of \f$ L_h u_h \f$
     */
    Precision apply(
        const DiscreteFunction& u,
        const Position pos,
        const Index sx,
        const Index sy) const;
        
    /**
     * \brief returns the coefficient of the center element at the given point
     * 
     * getCenter returns the center element of the stencil. Because for some
     * stencils it can make a diffrence if they are applied to a point on the
     * boarder or the center there is the parameter Position.
     * \see Position
     * \see apply
     * 
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position 
     * \param[in] sy        the y position
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the center coefficient
     */
    Precision getCenter(
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Precision hx,
        const Precision hy,
        const Point origin) const;
        
    /**
     * \brief retruns the coefficents of the stencil at the given point
     * 
     * getL() returns the coefficients of the stencil as valarray at the
     * given point. L[0] is always the coefficient of the center element
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
     * (This matrix shows the numbering of all getL return values.) 
     * All other elements are stored in unspecific order. To get the coordinates
     * of the other entries use getJx and getJy.\n
     * A Stencil also have to give the garanty that the elements 0-4 are all set
     * if size is one and that all values mentioned in the above matrix are set
     * if the size is two.\n
     * Because for some stencils it can make a diffrence if they are applied
     * to a point on the boarder or the center set the paramter Position correctly.\n
     * E.g. Laplacian2D2 with the stepsize hx=hy=1 and \f$a_x = a_y =1.0\f$
     * looks like:
     * \f[
     * \left[\begin{array}{ccc}
     *      &-1 &   \\
     * -1   & 4 &   -1\\
     *      &-1 &   \\
     * \end{array}\right]
     * \f]
     * So we have:\n
     * L    = {4,-1,-1,-1,-1}\n
     * Jx   = {0,-1,0,1,0}\n
     * Jy   = {0,0,1,0,-1}\n
     * \see getJx
     * \see getJy
     * \see Position
     * \see apply
     * \see Laplacian2D2
     * 
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position
     * \param[in] sy        the y position
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the coefficients of the stencil
     */
    virtual NumericArray getL(
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Precision hx,
        const Precision hy,
        const Point origin) const =0;
    
	/**
	 * \brief retruns the coefficents of the stencil at the given point
	 *
	 * 
	 *
	 * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position
     * \param[in] sy        the y position
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
	 * \param[in] size		the desired size of the returned stencil
     * \return              the coefficients of stencil in the desired size
     */
    NumericArray getLInSize(
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Precision hx,
        const Precision hy,
        const Point origin,
		const Index size ) const;

    /**
     * \brief returns the coordinate vector in x dir. for the coefficient vector
     * \see getL    for a more detailed description
     * 
     * \param[in] pos  the relative position to get the coordinates vector at
     * \param[in] nx   number of points in x direction (stride)
     * \param[in] ny   number of points in y direction
     * \return         the coordinates vector in x direction
     */
    virtual PositionArray getJx(
        const Position pos,
        const Index nx,
        const Index ny ) const =0;
    
    /**
     * \brief returns the coordinate vector in y dir. for the coefficient vector
     * \see getL    for a more detailed description
     * 
     * \param[in] pos  the relative position to get the coordinates vector at
     * \param[in] nx   number of points in x direction (stride)
     * \param[in] ny   number of points in y direction
     * \return          the coordinates vector in y direction
     */
    virtual PositionArray getJy(
        const Position pos,
        const Index nx,
        const Index ny ) const =0;
    
    /**
     * \brief pushs new TransferOperators on the stack
     * 
     * pushTransferOperators is needed for galerkin like operators. Each time
     * we are going to a coarser grid we need to tell the galerkin operator
     * what TransferOperators we are using.
     * 
     * \param[in] restriction   the Restriction to put on the stack 
     * \param[in] prolongation  the Prolongation to put on the stack
     * \param[in] nx            the number of steps in x-dir. on coarse grid
     * \param[in] ny            the number of steps in y-dir. on coarse grid
     */
    virtual void pushTransferOperators(
        const Restriction&,
        const Prolongation&,
        const Index nx,
        const Index ny);
    
    /**
     * \brief removes the top TransferOperators from the stack
     * 
     * \see pushTransferOperators
     */
    virtual void popTransferOperators();
    
    /**
     * \brief returns the size of the stencil
     * 
     * The size of the stencil is the max entry of the valarrays returnd by
     * getJx and getJy.
     * \see getJx
     * \see getJy
     * 
     * \return  the size of the stencil
     */
    virtual Index size() const =0;
    
    /**
     * \brief isConstant says if the stencil has constant coefficients
     * 
     * \return      true if the stencil has constant coefficeint
     */
    virtual bool isConstant() const =0;
private:
	//we don't want the autogenerated copy constructor and assignment operator
	Stencil( const Stencil& );
    Stencil& operator= ( const Stencil& );
};

}

#endif /*STENCIL_H_*/
