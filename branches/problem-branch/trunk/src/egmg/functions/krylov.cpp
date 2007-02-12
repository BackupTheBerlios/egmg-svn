/*
 *  krylov.cpp
 *  
 *
 *  Created by Matthias Rettenmeier on 27.08.06.
 *  Copyright 2006. All rights reserved.
 *
 */
#include <ostream>
#include <iomanip>
#include <iostream>
#include <valarray>
#include "krylov.h"

namespace mg
{
	NumericArray krylov(
		const size_t m,
		const NumericArray& Defect,
		const Index nx,
		const Index ny)
	{
		if (Defect.size()!=m*(nx+1)*(ny+1))
			{
			std::cout<<"Krylov can't proceed due to incorrect amount of Data"<<std::endl;
			}
		if (m<2)
			{
			std::cout<<"To little information! Can't proceed"<<std::endl;
			}
		//-- Setup linear System with Krylovspaces ( Ax=Vec ) 
		NumericArray Vec(0.0, m-1);
		NumericArray A=createKS(m, Defect, Vec, nx, ny);
		//-- Use QR-Dismantling to solve System ( QRx=Vec )
		NumericArray Q(0.0, ((m-1)*(m-1)));
		for(size_t i=0; i<(m-1); i++) Q[i*m]=1;
		DismantlingQR(m-1, 0, A, Q);
		//-- Multiply right side with Q^T and solve linear system ( Rx=Q^T * Vec )		
		//-- here: Q^T = Q
		NumericArray Vec2 = productVec(Q, Vec, m-1);
		NumericArray Koef = solveSys(A, Vec2, m-1);
		return Koef;
	}

	NumericArray createKS(
		int m,
		const NumericArray& Defect,
		NumericArray& Vec,
		const Index nx,
		const Index ny)
	{	
		const Index size=(nx+1)*(ny+1);
		NumericArray result(0.0,(m-1)*(m-1));
		NumericArray Vec1(0.0,size);
		size_t points=(nx-1)*(ny-1);

		//-- get defect vectors from defect history
		//-- and setup KrylovSystem 
		std::slice V1((m-1)*size,size,1);
		Vec1=Defect[V1];
		for (int i=0; i<(m-1); i++)
		{
			NumericArray Vec2(0.0, size);
			std::slice V2((m-2-i)*size,size,1);
			Vec2=Defect[V2];
			for(int k=0; k<(m-1); k++)
			{
				NumericArray Vec3(0.0, size);
				std::slice V3((m-2-k)*size,size,1);
				Vec3=Defect[V3];
				result[i*(m-1)+k]
				  =(Vec2 * Vec3).sum()/points
					-(Vec1 * Vec2).sum()/points
					-(Vec1 * Vec3).sum()/points
					+(Vec1 * Vec1).sum()/points;
			}
		}
		for(int i=0; i<(m-1); i++)
		{
			NumericArray Vec4(0.0, size);
			std::slice V4((m-2-i)*size,size,1);
			Vec4=Defect[V4];
			Vec[i]
				=(Vec1*Vec1).sum()/points
				-(Vec1*Vec4).sum()/points;
		}
		return result;
	}

	void DismantlingQR(
		size_t dim,
		size_t depth,
		NumericArray& A,
		NumericArray& T)
	{
		size_t start=(depth*dim)+depth;
		size_t newdim=dim-depth;
		Precision sign=1;
		NumericArray H(0.0, (dim*dim));
		NumericArray h(0.0, (newdim));
		std::valarray<size_t> sizes(2);
		std::valarray<size_t> strides(2);
		sizes[0]=newdim;
		sizes[1]=newdim;
		strides[0]=dim;
		strides[1]=1;
		std::gslice sub(start,sizes,strides);
		std::slice col1(start,newdim,dim);
			
		// create "Vector" from first column
		h=A[col1];
		if(h[0]<0) sign=-1;
		h[0]+=sign*sqrt((h*h).sum());
		h*=(1/sqrt((h*h).sum()));
		// start of with H as Identity I
		for(size_t i=0; i<dim; i++) H[i*(dim+1)]=1.0;
		// look at sub_matrix from H
		NumericArray Temp = H[sub];
		if(newdim>1)
		{
			for(size_t i=0; i<newdim; i++) Temp[i*(newdim+1)]=h[i]*h[i];
		}
		else Temp[0]=0;
		for(size_t i=1; i<newdim; i++)
		{
			for(size_t j=0; j<(newdim-i); j++)
			{
				Temp[(newdim)*i+(newdim+1)*j]=Temp[i+(newdim+1)*j] = h[j]*h[i+j]; 
			}	
		}
		Temp*=2;
		H[sub]-=Temp;
		
		// calculate new A
		A=product(H, A, dim);
		T=product(H, T, dim);
		if (newdim>1) DismantlingQR(dim, depth+1, A, T);
	}
	
	NumericArray product(
		NumericArray& A, 
		NumericArray& B,
		size_t dim)
	{
		NumericArray result(0.0, dim*dim);
		NumericArray Arow(0.0, dim);
		NumericArray Bcol(0.0, dim);
		for(size_t i=0; i< dim; i++)
		{
			std::slice ROW(i*dim, dim, 1);
			Arow=A[ROW];
			for(size_t j=0; j<dim; j++)
			{
				std::slice COLUMN(j, dim, dim);
				Bcol=B[COLUMN];
				result[i*dim+j]=(Arow*Bcol).sum();
			}
		}
		return result;
	}
	
	NumericArray productVec(
		NumericArray& A, 
		NumericArray& B,
		size_t dim)
	{
		NumericArray result(0.0, dim);
		NumericArray Arow(0.0, dim);
		for(size_t i=0; i< dim; i++)
		{
			std::slice ROW(i*dim, dim, 1);
			Arow=A[ROW];
			result[i]=(Arow*B).sum();
		}
		return result;
	}
	
	NumericArray solveSys(
		NumericArray A, 
		NumericArray Vec, 
		size_t dim)
	{
		//-- backward solving linear system
		NumericArray koef(0.0, dim);
		koef[dim-1]=Vec[dim-1]/A[dim*dim-1];
		for(int i=dim-2; i>=0; i--)
		{
			koef[i]=Vec[i];
			for(int j=dim-1; j>i; j--)
			{
				koef[i]-=koef[j]*A[(i*dim)+j];
			}
			koef[i]=koef[i]/A[(i*dim)+i];
		}
		return koef;
	}
		
	void KrylovUpdate(
		NumericArray& u,
		NumericArray& UHistory,
		NumericArray& Koef,
		size_t size,
		size_t dim)
	{
		//-- adding linear combination to most recent multigrid solution
		std::slice last(dim*size, size, 1);
		NumericArray Um=UHistory[last];
		NumericArray Temp(0.0, size);
		for(size_t  i=0; i< dim; i++)
		{
			std::slice KR( (dim-1-i)*size,size,1 );
			Temp=UHistory[KR];
			u+=Koef[i]*(Temp-Um);
		}
	}
	
		
}

