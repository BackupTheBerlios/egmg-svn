/** \file parameters.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief parameters.h provides typedefs used in the whole program
 */
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include<valarray>

/**
 * \brief the namespace mg is for everything in this library
 * 
 * mg is an abbreviation for multigrid
 */
namespace mg
{
	/**
	 * \brief the precision to do the calculations with
	 */
	typedef double precision;
	
	/**
	 * \enum mg::pos
	 * \brief pos describes relative positions
	 * 
	 * The meaning of the constants are:
	 * <table border="0">
	 * <tr>
	 * <td>c</td><td>center</td>
	 * </tr><tr>
	 * <td>w</td><td>west</td>
	 * </tr><tr>
	 * <td>nw</td><td>north west</td>
	 * </tr><tr>
	 * <td>n</td><td>north</td>
	 * </tr><tr>
	 * <td>ne</td><td>north east</td>
	 * </tr><tr>
	 * <td>e</td><td>east</td>
	 * </tr><tr>
	 * <td>se</td><td>south east</td>
	 * </tr><tr>
	 * <td>s</td><td>south</td>
	 * </tr><tr>
	 * <td>sw</td><td>south west</td>
	 * </table>
	 */
	enum pos{c=0,w=1,n=2,e=3,s=4,nw=5,ne=6,se=7,sw=8};
	
	/**
	 * \enum mg::direction
	 * \brief direction describes direction of line relaxation
	 * 
	 * The meaning of the constants are:
	 * <table border="0">
	 * <tr>
	 * <td>alt_zebra</td><td>alternating zebra</td>
	 * </tr><tr>
	 * <td>xline</td><td>x direction</td>
	 * </tr><tr>
	 * <td>yline</td><td>y direction</td>
	 * </tr><tr>
	 * <td>x_zebra</td><td>zebra x direction</td>
	 * </tr><tr>
	 * <td>y_zebra</td><td>zebra y direction</td>
	 * </tr><tr>
	 * </table>
	 */
	enum direction{altline = 0, xline = 1, yline =2, alt_zebra=3, x_zebra = 4, y_zebra = 5};
	
	/**
	 * \brief a function pointer to a function \f$ f: K^2 \rightarrow K \f$
	 */
	typedef precision (*function2D) (const precision x, const precision y);
	
}
#endif /*PARAMETERS_H_*/
