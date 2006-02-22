/** \file output.cpp
 * \author Matthias Rettenmeier
 * \brief contains the implementation of the functions output and outputvalue
 */
#include "output.h"
#include<iomanip>
namespace mg
{

int output(std::vector<Precision>& vec, size_t flag){
  Precision temp0; //temp variable for residnull
  Precision tempfive; //temp variable for residuen after 5th cycle.
  Precision temp1, temp2;
  size_t max=(size_t)vec.size(); //values stored in history

  //error check
  if (vec.size() <= 0){
    std::cout << "Error accured! Vector empty!" << std::endl;
    return -1;
  }

  // create three arrays to store the calculated convergence rates
  std::valarray<Precision> F1(0.0,max);
  std::valarray<Precision> F2(0.0,max);
  std::valarray<Precision> F3(0.0,max);

  // initializing temp0, tempfive
  temp0 = vec.at(0);
  if (vec.size() >=5 ) tempfive = vec.at(5);
  else tempfive = -1;

  // decision on headline
  if(flag==1) std::cout << "Defektreduktionsfaktoren bzgl. der Max-Norm: "<< std::endl;
  else if(flag==2) std::cout << "Defektreduktionsfaktoren bzgl. der 2-Norm: " <<std::endl;
  else std::cout << "Defektkonvergenzfaktoren (Norm nicht definiert): " <<std::endl;

  // go through the history of residues and calculate convergence rates for each cycle
  std::cout<<"F1=(||r_i||/||r_o||)^(1/i)\nF2=(||r_i||/||r_5||)^(1/(i-5))\nF3=||r_i||/||r_(i-1)||"<<std::endl;
  for(size_t i=1; i<max; i++){
    temp1 = vec.at(i);
    temp2 = vec.at(i-1);
    
    F1[i]= pow( (temp1/temp0) , static_cast<double> (1./i));
	//change by jiri old:
 	//   if(i>5 && tempfive != -1 ) F2[i]=pow( (temp1/tempfive) , static_cast<Precision> (1./(i-5)) );
	//new:
	if(i>5 && vec.size() >= 5 ) F2[i]=pow( (temp1/tempfive) , static_cast<Precision> (1./(i-5)) );
    else F2[i]=-1;

    F3[i]= temp1/temp2;
    
    // output the calculated content
    std::cout<<"F1 = "<<std::setw(10)<<F1[i] <<" | F2 = "<<std::setw(10)<< F2[i] << " | F3 = "<<std::setw(10)<< F3[i]<<std::endl;
  }
  return 0;
}



int outputvalue(std::valarray<Precision> u, size_t Nx, size_t Ny){
  size_t x,y;

   //error check
  if (u.size() <= 0 || u.size()!= (Nx+1)*(Ny+1)){
    std::cout << "Error accured! Output aborted." << std::endl;
    return -1;
  }

  std::cout<<"Welcher Punkt soll ausgegeben werden?"<<std::endl;
  std::cout<<"X: " ;
  std::cin >> x;
  std::cout<<"Y: " ;
  std::cin >> y;

  if( x > Nx || y > Ny ){
    std::cout<<"Wert liegt außerhalb des Gitters"<<std::endl;
    outputvalue(u, Nx, Ny);
  }
  else std::cout<<"Ausgabe für Punkt("<< x <<","<< y <<"): "<< u[x+(Nx+1)*y] <<std::endl;
  
  return 0;
}

}
