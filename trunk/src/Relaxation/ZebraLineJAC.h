/** \file lineGS.h
 * \author Andre Oeckerath
 * \see Relaxation.h
 */
#ifndef lineJAC_H_
#define lineJAC_H_

#include "Relaxation.h"
#include "../Stencil/Stencil.h"


namespace mg
{
/**
 * \brief lineJAC is a class template for a Jacobi line relaxation
 */

	class ZebraLineJAC : public mg::Relaxation
	{
        private:
		
	        int nreld_;
	        int nrelu_;
			Precision omega_;
			int direction_;

			void ninepointxline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void ninepointyline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void ninepointxzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void ninepointyzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void xline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void yline(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void xzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

			void yzebra(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
		                    std::valarray<Precision> resid, const Stencil &stencil, const size_t Nx, 
					        const size_t Ny) const;

		        

	    public:
            /**
	         * \brief The constructor of a lineJAC object
	         * 
	         * lineJAC constructs lineJAC object with:
	         * \param nreld		number of pre smothing steps (default 1)
	         * \param nrelu		number of post smothing steps (default 1)
			 * \param direction_ 
             * altline = alternating x/y (default 0)
			 * xline = x line
             * yline = y line
			 * alt_zebra = alternating x/y zebra
			 * x_zebra = x line zebra
			 * y_zebra = y line zebra
             */ 
	        ZebraLineJAC(const int nreld =1,const int nrelu =1, const Precision omega =1.0,const int direction =0) : nreld_(nreld), nrelu_(nrelu), omega_(omega), direction_(direction) {}
            virtual ~ZebraLineJAC() {}
		    
			/**
	          * \brief getPreSmoothingSteps() returns the number of pre smoothing steps
	          * \return the number of pre somthing steps
	          */
		    int get_nreld() const {return nreld_;}
            void set_nreld(int nreld) {nreld_=nreld;}

			/**
	          * \brief getPostSmoothingSteps() returns the number of post smoothing steps
	          * \return the number of post somthing steps
	          */
		    int get_nrelu() const {return nrelu_;}
	        void set_nrelu(int nrelu) {nrelu_=nrelu;}

			
            /**
	          * \brief get_direction() returns direction of line smoothing
	          * 0 = alternating x/y
			  * 1 = x line
              * 2 = y line
			  * 3 = alternating x/y zebra
			  * 4 = x line zebra
			  * 5 = y line zebra
	          */
			int get_direction() const {return direction_;}
			void set_direction(int direction) {direction_=direction;}
			/**
	          * \brief relax() executes one relaxation step on the input vector
	          * 
	          * relax() exectues one Jacobi line relaxation step on the 
			  * input vector on a rectangular 2D gird with lexicographic ordering and the
	          * discretazation Stencil with const coefficient for a pde
	          * \param u		the vector representation of the 2D grid to perform the
	          * 				relaxation on this vector will be changed
	          * \param fv		the right hand side of the pde
	          * \param Nx	number of steps in x direction
	          * \param Ny	number of steps in y direction
			  */
            void relax(std::valarray<Precision> &u, const std::valarray<Precision> &fv, 
				const Stencil &stencil, const size_t Nx, const size_t Ny) const;
		
    };

}



#endif
