/*
 *  krylov.h
 *  
 *
 *  Created by Matthias Rettenmeier on 27.08.06.
 *
 */
#ifndef KRYLOV_H_
#define KRYLOV_H_

#include "../general/parameters.h"
#include <algorithm>

namespace mg
{
/**
 * \brief krylov() calculates coefficients needed for krylov accaleration
 *
 * krylov() requires a history of defects to setup a linear system which 
 * is then solved to receive coeffictients for the krylov acceleration. 
 * 
 * \param m				Size of the krylov space that shall be used
 * \param[in] Defects	NumericArray that includes a history of defects
 *						required to setup Krylov space. The number of defects
 *						must equal m.
 * \param[in] nx        number of steps in x direction
 * \param[in] ny        number of steps in y direction
 * \return				(m-1) coefficients for krylov accaleration
 */
NumericArray krylov(
	const size_t m,
	const NumericArray& Defects,
	const Index nx,
	const Index ny);
/**
 * \brief createKS() sets up the linear system which is solved to receive krylov coefficients
 * 
 * \param[in] m				Size of the krylov space that shall be used
 * \param[in, out] Defects	NumericArray that includes a history of defects
 *							required to setup Krylov space. The number of defects
 *							must equal m.
 * \param[in, out] Vec		NumericArray that contains right side of linear system
 * \param[in] nx			number of steps in x direction
 * \param[in] ny			number of steps in y direction
 * \return					NumericArray that contains left side of linear system 
 */
NumericArray createKS(
	int m,
	const NumericArray& Defects,
	NumericArray& Vec,
	const Index nx,
	const Index ny);
/**
 * \brief DismatlingQR() dismatles given matrix A into an orthogonal and an upper triagle matrix
 * 
 * \param[in] dim			Dimension of matrix A
 * \param[in, out] depth	Internal variable used to control recursive calls
 * \param[in, out] A		Matrix which shall be dismatled. 
 *							Will be transformed to upper triangle matrix R
 * \param[in, out] T		Unit matrix that will be transformed to orthogonal matrix Q^T
 */
void DismantlingQR(
		size_t dim,
		size_t depth,
		NumericArray& A,
		NumericArray& T);
/**
 * \brief product() Function to calculate the product between two (nxn)-matrixes
 * 
 * \param[in] A		(n*n)-matrix
 * \param[in] B		(n*n)-matrix
 * \param[in] dim	n
 * \return			result of multiplication
 */	
NumericArray product(
		NumericArray& A, 
		NumericArray& B,
		size_t dim);
/**
 * \brief productVec() Function to calculate the product between a (nxn)-matrix and a (n)-vector
 * 
 * \param[in] A		(n*n)-matrix
 * \param[in] B		(n)-vector
 * \param[in] dim	n
 * \return			result of multiplication
 */
NumericArray productVec(
		NumericArray& A, 
		NumericArray& B,
		size_t dim);
/**
 * \brief solveSys() solves linear system backward
 * 
 * \param[in] A		upper triagle matrix (left side of linear system)
 * \param[in] Vec	vector (right side of linear system)
 * \param[in] dim	dimension of linear system
 * 
 * \return			solution to linear system
 */
NumericArray solveSys(
		NumericArray A, 
		NumericArray Vec, 
		size_t dim);

/**
 * \brief KrylovUpdate()	applies linear combination of previous aquired solutions 
 *							to most recent solution in respect to krylov optimization.
 * 
 * \param[in, out] u		most recent solution of multigrid itteration
 * \param[in] UHistory		history of multigrid solutions 
 * \param[in] Koef			previously calculated Krylov coefficients
 * \param[in] size			(nx+1)*(ny+1)
 * \param[in] dim			size of Krylov space previously used (m)
 */

void KrylovUpdate(
		NumericArray& u,
		NumericArray& UHistory,
		NumericArray& Koef,
		size_t size,
		size_t dim);


}
#endif /*KRYLOV_H_*/

