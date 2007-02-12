/** \file StencilBoard.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface the class StencilBoard.
 */
#ifndef STENCILBOARD_H_
#define STENCILBOARD_H_

#include "../general/parameters.h"
#include <vector>
#include <utility>
#include <map>


namespace mg
{

/**
 * \brief StencilBoard is a container for coefficients of differential operators
 *        in stencil notation.
 * 
 * StencilBoard stores the coefficients and position arrays of stencils acording
 * to the Stencil interface (\see Stencil). It is used e.g. by the class
 * Galerkin (\see Galerkin).
 */
class StencilBoard
{
public:
    /**
     * \brief Constructs an empty StencilBoard.
     */
	StencilBoard();

	virtual ~StencilBoard();
    
    /**
     * \brief insert a stencil into StencilBoard.
     * 
     * insert a stencil into the StencilBoard. The coefficient array
     * is a NumericArray with the coefficients of the stencil. Together with two
     * position arrays it describes the stencil. (\see Stencil.getL)
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position of the center element 
     * \param[in] sy        the y position of the center element
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \param[in] opL       the coefficient vector of the 
     */
    void insert( const Index level,
                 const Position pos,
                 const Index sx,
                 const Index sy,
                 const NumericArray& opL );
                 
    /**
     * \brief insert a stencil into StencilBoard.
     * 
     * insert a stencil into the StencilBoard. The position arrays describe,
     * together with the coefficient array, the stencil. (\see Stencil.getL)
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \param[in] jX        the position array in x direction
     * \param[in] jY        the position array in y direction
     */
    void insert( const Index level,
                 const Position pos,
                 const PositionArray& jX,
                 const PositionArray& jY );
                 
    /**
     * \brief check if this stencil board contains a stencil
     * 
     * contains checks if this stencil board contains the stencil at
     * at:
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position of the center element 
     * \param[in] sy        the y position of the center element
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              true if the coefficient vector was previously
     *                      inserted
     */
    bool contains( const Index level,
                   const Position pos,
                   const Index sx,
                   const Index sy ) const;
                   
    /**
     * \brief gets the coefficient vector at the given position
     * 
     * getL returns the coefficient vector at the given position.\n
     * CAVEATE: getL does not do a range check use contains for this if
     * necessary.
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] sx        the x position of the center element 
     * \param[in] sy        the y position of the center element
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the coefficient vector
     */
    const NumericArray& getL( const Index level,
                              const Position pos,
                              const Index sx,
                              const Index sy ) const;
    /**
     * \brief gets the coefficient vector in x dir. at the given position
     * 
     * getJx returns the position array in x direction at the given position.\n
     * CAVEATE: getJx does not do a range check use contains for this if
     * necessary.
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the position array in x dir
     */
    const PositionArray& getJx( const Index level,
                                const Position pos ) const;
                            
    /**
     * \brief gets the coefficient vector in y dir. at the given position
     * 
     * getJy returns the position array in y direction at the given position.\n
     * CAVEATE: getJy does not do a range check use contains for this if
     * necessary.
     * \param[in] pos       relative position in the domain according to the
     *                      enum Postion
     * \param[in] nx        number of points in x direction (stride)
     * \param[in] ny        number of points in y direction
     * \return              the position array in y dir
     */
    const PositionArray& getJy( const Index level,
                                const Position pos ) const;

private:
    //we don't want the autogenerated copy constructor and assignment operator
    StencilBoard( const StencilBoard& );
    StencilBoard& operator=( const StencilBoard& );
    
    typedef std::map<std::pair<Index,Index>,NumericArray> CoefficientMap;
    
    std::vector<CoefficientMap >                        m_coefficients[9];
    std::vector<PositionArray>                          m_posJx[9];
    std::vector<PositionArray>                          m_posJy[9];
};

}

#endif /*STENCILBOARD_H_*/
