/** \file parameters.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief parameters.h provides typedefs used in the whole program
 */
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

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
	 * \enum mg::Postion
	 * \brief Postion describes relative positions
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
	enum Position{c=0,w=1,n=2,e=3,s=4,nw=5,ne=6,se=7,sw=8};
	
	/**
	 * \brief a function pointer to a function \f$ f: K^2 \rightarrow K \f$
	 */
	typedef Precision (*function2D) (
		const Precision x,
		const Precision y);
	
}
#endif /*PARAMETERS_H_*/
