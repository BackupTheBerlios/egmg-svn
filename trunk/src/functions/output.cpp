/** \file output.cpp
 * \author Matthias Rettenmeier
 * \brief contains the implementation of the functions output and outputvalue
 */
#include <iomanip>
#include <cmath>
#include <iostream>
#include <ostream>
#include <stdexcept>

#include "output.h"

namespace mg
{

void convergenceRates(const std::vector<Precision>& vec, std::ostream& out)
{
    Precision resNull=vec.at(0);
    Precision resFive=-1; //temp variable for residuen after 5th cycle
    if (vec.size() >=5 )
        resFive=vec.at(5);

    // go through the history of residues and calculate convergence rates for 
    // each cycle
    out<<"CFI = (||r_i||/||r_o||)^(1/i)"<<std::endl
       <<"CF5 = (||r_i||/||r_5||)^(1/(i-5))"<<std::endl
       <<"CF  = ||r_i||/||r_(i-1)||"<<std::endl;
    out<<std::setw(10)<<"CFI"<<" | "
       <<std::setw(10)<<"CF5"<<" | "
       <<std::setw(10)<<"CF"<<std::endl;
    for(Index i=1; i<vec.size(); i++)
    {
        Precision tempRes1=vec.at(i);
        Precision tempRes2=vec.at(i-1);
        Precision residFracI=pow(tempRes1/resNull,static_cast<double>(1./i));
        Precision residFracV=-1;
        if(i>5 && vec.size() >= 5 )
            residFracV=pow(tempRes1/resFive,static_cast<Precision>(1./(i-5)));
        Precision residFrac=tempRes1/tempRes2;
        // output the calculated content
        out<<std::setw(10)<<residFracI<<" | "
           <<std::setw(10)<<residFracV<<" | "
           <<std::setw(10)<<residFrac<<std::endl;
    }
}

void gnuPlotDiscreteFunction(
    const NumericArray& u,
    Index nx,
    Index ny,
    std::ostream& out)
{
    if ((nx+1)*(ny+1)!=u.size())
        throw std::domain_error("u");
    out<<"#begin points"<<std::endl;
    out<<std::setw(10)<<std::left<<"#x"<<" "
       <<std::setw(10)<<std::left<<"y"<<" "
       <<std::setw(10)<<std::left<<"value"<<std::endl;
    Precision hx=1.0/nx;
    Precision hy=1.0/ny;
    for (Index sy=0; sy<=ny; ++sy)
        for (Index sx=0; sx<=nx; ++sx)
            out<<std::setw(10)<<std::left<<sx*hx<<" "
               <<std::setw(10)<<std::left<<sy*hy<<" "
               <<std::setw(10)<<std::left<<u[sy*(nx+1)+sx]<<std::endl;
    out<<"#end points"<<std::endl;
}
}
