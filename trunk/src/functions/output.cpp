/** \file output.cpp
 * \author Matthias Rettenmeier
 * \brief contains the implementation of the functions output and outputvalue
 */
#include<iomanip>
#include<cmath>
#include<iostream>

#include "output.h"

namespace mg
{

void convergenceRates(std::vector<Precision>& vec)
{
    Precision resNull=vec.at(0);
    Precision resFive=-1; //temp variable for residuen after 5th cycle
    if (vec.size() >=5 )
        resFive=vec.at(5);

    // go through the history of residues and calculate convergence rates for 
    // each cycle
    std::cout<<"CFI = (||r_i||/||r_o||)^(1/i)"<<std::endl
             <<"CF5 = (||r_i||/||r_5||)^(1/(i-5))"<<std::endl
             <<"CF  = ||r_i||/||r_(i-1)||"<<std::endl;
    std::cout<<std::setw(10)<<"CFI"<<" | "
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
        std::cout<<std::setw(10)<<residFracI<<" | "
                 <<std::setw(10)<<residFracV<<" | "
                 <<std::setw(10)<<residFrac<<std::endl;
    }
}
}
