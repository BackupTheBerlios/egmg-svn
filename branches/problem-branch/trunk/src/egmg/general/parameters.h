/** \file parameters.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief parameters.h provides typedefs used in the whole program
 */
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <valarray>
#include <cassert>
#include <cmath>
#include <utility>
#include <memory>

/**
 * \brief the namespace mg is for everything in this library
 * 
 * mg is an abbreviation for multigrid
 */
namespace mg
{
    class Problem;
    
    typedef std::auto_ptr<Problem> ProblemPtr;
    /**
     * \brief the Precision to do the calculations with
     */
    typedef double Precision;

    const Precision Pi=4*std::atan(static_cast<Precision>(1.0));
    
    /**
     * \brief the index
     */
    typedef size_t Index;

    /**
     * \brief Integer Values
     */
    typedef int Integer;
    
    /**
     * \brief vector holding numeric values
     */
    typedef std::valarray<Precision> NumericArray;
    
    /**
     * \brief vector holding integers
     */
    typedef std::valarray<Integer> PositionArray;
    
    /**
     * \brief pair of matrix and right and side
     */
    typedef std::pair<NumericArray,NumericArray> LinearEquationSystem;
        
    /**
     * \enum mg::Position
     * \brief Position describes relative positions
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
     * </tr>
     * </table>
     */
    enum Position{C=0,W=1,N=2,E=3,S=4,NW=5,NE=6,SE=7,SW=8};
                 

    enum PositionExtension {WW=9,NWW=10,NNWW=11,NNW=12,
        NN=13,NNE=14,NNEE=15,NEE=16,EE=17,SEE=18,
        SSEE=19,SSE=20,SS=21,SSW=22,SSWW=23,SWW=24,NamedPositions=25};
    
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
     * </tr>
     * </table>
     */
    enum Direction{ALTDIR=0,XDIR=1,YDIR=2};
    
    struct IndexPair
    {
        IndexPair(Index sx, Index sy): sx(sx), sy(sy){};
        Index sx;
        Index sy;
    };
    
    struct Point
    {
        Point(Precision x, Precision y): x(x), y(y) {};
        Precision x;
        Precision y;
    };
    
    bool operator== (const Point& lhs, const Point& rhs );

    #ifndef NASSERT
        #define ASSERT(exp) assert(exp);
    #else
        #define ASSERT(exp) ((void)0);
    #endif
    
}

#endif /*PARAMETERS_H_*/
