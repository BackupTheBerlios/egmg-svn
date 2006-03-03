/** \file parameters.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief parameters.h provides typedefs used in the whole program
 */
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <valarray>

/**
 * \brief the namespace mg is for everything in this library
 * 
 * mg is an abbreviation for multigrid
 */
namespace mg
{
	/**
	 * \brief the Precision to do the calculations with
	 */
	typedef double Precision;
    
    /**
     * \brief the index
     */
    typedef size_t Index;
    
    /**
     * \brief vector holding numeric values
     */
    typedef std::valarray<Precision> NumericArray;
    
    /**
     * \brief vector holding integers
     */
    typedef std::valarray<int> PositionArray;
		
	/**
	 * \enum mg::Postion
	 * \brief Postion describes relative positions
	 * 
	 * The meaning of the constants are:
	 * <table border="0">
	 * <tr>
	 * <td>C</td><td>center</td>
	 * </tr><tr>
	 * <td>W</td><td>west</td>
	 * </tr><tr>
	 * <td>NW</td><td>north west</td>
	 * </tr><tr>
	 * <td>N</td><td>north</td>
	 * </tr><tr>
	 * <td>NE</td><td>north east</td>
	 * </tr><tr>
	 * <td>E</td><td>east</td>
	 * </tr><tr>
	 * <td>SE</td><td>south east</td>
	 * </tr><tr>
	 * <td>S</td><td>south</td>
	 * </tr><tr>
	 * <td>SW</td><td>south west</td>
	 * </table>
	 */
	enum Position{C=0,W=1,N=2,E=3,S=4,NW=5,NE=6,SE=7,SW=8};
    
    /**
     * \enum mg::Direction
     * \brief Direction describes direction of line relaxation
     * 
     * The meaning of the constants are:
     * <table border="0">
     * <tr>
     * <td>ALTDIR</td><td>alternating directions</td>
     * </tr><tr>
     * <td>XDIR</td><td>x direction</td>
     * </tr><tr>
     * <td>YDIR</td><td>y direction</td>
     * </tr><tr>
     * </table>
     */
    enum Direction{ALTDIR=0,XDIR=1,YDIR=2};
	
	/**
	 * \brief a function pointer to a function \f$ f: K^2 \rightarrow K \f$
	 */
	typedef Precision (*function2D) (
		const Precision x,
		const Precision y);
	
}
#endif /*PARAMETERS_H_*/
