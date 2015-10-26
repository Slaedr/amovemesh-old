/* Implemetation of sparse storage schemes for matrices.
   Aditya Kashi
   21 July 2015
*/

#define __ASPARSEMATRIX_H

// for memcpy() etc
#ifndef _GLIBCXX_CSTRING
#include <cstring>
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLBCXX_CMATH
#include <cmath>
#endif

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif

#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#define nthreads_sm 8
#endif
#endif

namespace amat {

typedef Matrix<double> Mat;

class SparseMatrix
/* Abstract class for sparse storage of matrices */
{
protected:
	double* val;
	int nnz;
	int nrows;
	int ncols;
	int valsize;
public:
	SparseMatrix() {}
	SparseMatrix(int num_rows, int num_cols) : nrows(num_rows), ncols(num_cols)
	{
		nnz = 0;
	}
	int rows() { return nrows; }
	int cols() { return ncols; }
	virtual void set(int x, int y, double value) = 0;
	virtual double get(int x, int y) = 0;
	virtual void mprint() = 0;

	virtual void multiply(Mat& x, Mat* a, char paralel) = 0;
	virtual void multiply_parts(Mat* x, Mat* y, Mat* ans, int p) = 0;
};

class MatrixCOO : public SparseMatrix
/*
 Implements sparse matrix storage in coordinate format. val contains non-zero values, col_ind and row_ind contain corresponding
 column and row indices; all three have length nnz after all values are stored.
*/
{
	int* col_ind;
	int* row_ind;
	bool sorted;		// a boolean flag that is true if row_ind is sorted in increasing order, and for each row index, col_ind is sorted in increasing order
	const double tol = 1e-60;			// tolerance below which to consider an entry zero
public:
	MatrixCOO() { val = new double; row_ind = new int; col_ind = new int; valsize = 1; sorted = false; }

	MatrixCOO(int num_rows, int num_cols) : SparseMatrix(num_rows, num_cols)
	{
		// Size of the 3 arrays is initially 2*num_rows

		valsize = 2*nrows;
		val = new double[valsize];
		col_ind = new int[valsize];
		row_ind = new int[valsize];
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0;
			col_ind[i] = -1;
			row_ind[i] = -1;
		}
		sorted = false;
		//tol = 1e-60;
	}

	MatrixCOO(const MatrixCOO& other)
	{
		valsize = other.valsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;
		sorted = other.sorted;
		val = new double[valsize];
		col_ind = new int[valsize];
		row_ind = new int[valsize];
		for(int i = 0; i < nnz; i++)
		{
			val[i] = other.val[i];
			row_ind[i] = other.row_ind[i];
			col_ind[i] = other.col_ind[i];
		}
		//tol = 1e-60;
	}

	MatrixCOO& operator=(MatrixCOO other)
	{
		valsize = other.valsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;
		sorted = other.sorted;

		delete [] val;
		delete [] row_ind;
		delete [] col_ind;

		val = new double[valsize];
		col_ind = new int[valsize];
		row_ind = new int[valsize];

		for(int i = 0; i < nnz; i++)
		{
			val[i] = other.val[i];
			row_ind[i] = other.row_ind[i];
			col_ind[i] = other.col_ind[i];
		}
		return *this;
	}

	void setup(int num_rows, int num_cols)
	{
		nrows = num_rows; ncols = num_cols;
		// Size of the 3 arrays is initially 2*num_rows
		valsize = 2*nrows;
		val = new double[valsize];
		col_ind = new int[valsize];
		row_ind = new int[valsize];
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0;
			col_ind[i] = -1;
			row_ind[i] = -1;
		}
		sorted = false;
		nnz = 0;
	}

	~MatrixCOO()
	{
		delete [] val;
		delete [] col_ind;
		delete [] row_ind;
	}

	void set(int x, int y, double value)
	{
		#pragma omp critical (omp_set)		// this block should be executed by only one thread at a time. 'omp_set' is the name of this critical block
		{
			if(row_ind[valsize-1] != -1)		// if last location of row_ind is not unused, allocate more space
			{
				std::cout << "MatrixCOO: set(): Increasing size of arrays\n";
				/*double* tval = new double[valsize];
				int* trow_ind = new int[valsize];
				int* tcol_ind = new int[valsize];

				for(int i = 0; i < valsize; i++)
				{
					tval[i] = val[i];
					trow_ind[i] = row_ind[i];
					tcol_ind[i] = col_ind[i];
				}

				delete [] val;
				delete [] row_ind;
				delete [] col_ind;

				val = new double[2*valsize];
				row_ind = new int[2*valsize];
				col_ind = new int[2*valsize];

				for(int i = 0; i < valsize; i++)
				{
					val[i] = tval[i];
					row_ind[i] = trow_ind[i];
					col_ind[i] = tcol_ind[i];
				}*/
				//std::cout << "Valsize is " << valsize << std::endl;
				double* tval = new double[2*valsize];
				int* trow_ind = new int[2*valsize];
				int* tcol_ind = new int[2*valsize];

				memcpy(tval, val, valsize*sizeof(double));
				memcpy(trow_ind, row_ind, valsize*sizeof(int));
				memcpy(tcol_ind, col_ind, valsize*sizeof(int));

				delete [] val;
				delete [] row_ind;
				delete [] col_ind;

				val = tval;
				row_ind = trow_ind;
				col_ind = tcol_ind;

				// set default values to newly created empty space
				for(int i = valsize; i < 2*valsize; i++)
				{
					val[i] = 0;
					row_ind[i] = -1;
					col_ind[i] = -1;
				}

				valsize *= 2;
			}

			// now add the new element if it's non-zero
			//std::cout << nnz << std::endl;
			if(value != 0.0)
			{
				int index = -1;
				for(int i = 0; i < nnz; i++)
				{
					//std::cout << i << " ";
					if(row_ind[i]==x)
					{
						if(col_ind[i]==y)
						{
							index = i;
							break;
						}
					}
				}
				if(index == -1)
				{
					val[nnz] = value;
					row_ind[nnz] = x;
					col_ind[nnz] = y;
					nnz++;
					sorted = false;
				}
				else
				{	// a non-zero element already exists in the position (x,y)
					val[index] = value;
				}
			}
		} // end omp critical
	}
	
	void set_fast(int x, int y, double value)
	{
		#pragma omp critical (omp_set)		// this block should be executed by only one thread at a time. 'omp_set' is the name of this critical block
		{
			if(row_ind[valsize-1] != -1)		// if last location of row_ind is not unused, allocate more space
			{
				std::cout << "MatrixCOO: set(): Increasing size of arrays\n";
				/*double* tval = new double[valsize];
				int* trow_ind = new int[valsize];
				int* tcol_ind = new int[valsize];

				for(int i = 0; i < valsize; i++)
				{
					tval[i] = val[i];
					trow_ind[i] = row_ind[i];
					tcol_ind[i] = col_ind[i];
				}

				delete [] val;
				delete [] row_ind;
				delete [] col_ind;

				val = new double[2*valsize];
				row_ind = new int[2*valsize];
				col_ind = new int[2*valsize];

				for(int i = 0; i < valsize; i++)
				{
					val[i] = tval[i];
					row_ind[i] = trow_ind[i];
					col_ind[i] = tcol_ind[i];
				}*/
				//std::cout << "Valsize is " << valsize << std::endl;
				double* tval = new double[2*valsize];
				int* trow_ind = new int[2*valsize];
				int* tcol_ind = new int[2*valsize];

				memcpy(tval, val, valsize*sizeof(double));
				memcpy(trow_ind, row_ind, valsize*sizeof(int));
				memcpy(tcol_ind, col_ind, valsize*sizeof(int));

				delete [] val;
				delete [] row_ind;
				delete [] col_ind;

				val = tval;
				row_ind = trow_ind;
				col_ind = tcol_ind;

				// set default values to newly created empty space
				for(int i = valsize; i < 2*valsize; i++)
				{
					val[i] = 0;
					row_ind[i] = -1;
					col_ind[i] = -1;
				}

				valsize *= 2;
			}

			// now add the new element if it's non-zero
			//std::cout << nnz << std::endl;
			if(value != 0.0)
			{
				val[nnz] = value;
				row_ind[nnz] = x;
				col_ind[nnz] = y;
				nnz++;
				sorted = false;
			}
		} // end omp critical
	}

	double get(int x, int y)
	{
		//std::cout << "MatrixCOO: get(): nnz is " << nnz << '\n';
		int index = -1;
		if(sorted == false)
		{
			int i;
			#ifdef _OPENMP
			int* row_ind = MatrixCOO::row_ind; int* col_ind = MatrixCOO::col_ind; int nnz = MatrixCOO::nnz;
			#endif
			#pragma omp parallel for default(none) private(i) shared(row_ind,col_ind,nnz,x,y,index) num_threads(nthreads_sm)
			for(i = 0; i < nnz; i++)
			{
				if(row_ind[i]==x)
				{
					if(col_ind[i]==y)
					{
						index = i;
						//break;
					}
				}
			}
			if(index == -1)
			{
				//std::cout << "! MatrixCOO: get(): Element at " << x << ", " << y << " not found!\n";
				return 0;
			}
			else return val[index];
		}
		else
		{
			// binary search
			int left = 0, right = nnz-1; int mid = (left+right)/2;
			while(left <= right)
			{
				mid = (left+right)/2;
				if(row_ind[mid] == x) break;
				else if(row_ind[mid] < x) left = mid+1;
				else right = mid-1;
			}
			//std::cout << "* get(): " << mid << '\n';

			if(row_ind[mid] != x)
				return 0;

			int i;
			bool inflag1 = false, inflag2 = false;		// flags to prevent infinite jumping between consecutive elements
			for(i = mid; i >= 0 && i < nnz && row_ind[i] == x;)
			{
				if(col_ind[i] == y)
				{
					index = i;
					break;
				}
				else if(col_ind[i] < y)
				{
					if(inflag2 == true) return 0;
					inflag1 = true;
					i++;
				}
				else
				{
					if(inflag1 == true) return 0;
					inflag2 = true;
					i--;
				}
			}
			if(index == -1)
			{
				//std::cout << "! MatrixCOO: get(): Element at " << x << ", " << y << " not found!\n";
				return 0;
			}
			return val[index];
		}
	}

	void row_sort(double* v, int* r, int* c, int len)
	{
		if(len==1) return;

		int len1, len2;
		len1 = len/2; len2 = len - len1;

		// NOTE: Can this be done in-place, perhaps using 1 or 2 temp locations per array?

		// But for now, copy the two halves of the arrays into newly-allocated arrays
		double* v1 = new double[len1];
		double* v2 = new double[len2];
		int* r1 = new int[len1];
		int* r2 = new int[len2];
		int* c1 = new int[len1];
		int* c2 = new int[len2];

		for(int i = 0; i < len1; i++)
		{
			v1[i] = v[i];
			r1[i] = r[i];
			c1[i] = c[i];
		}
		for(int i = len1; i < len; i++)
		{
			v2[i-len1] = v[i];
			r2[i-len1] = r[i];
			c2[i-len1] = c[i];
		}

		row_sort(v1,r1,c1,len1);
		row_sort(v2,r2,c2,len2);

		// now merge the two.

		int i = 0, j = 0, k = 0;
		while(k < len)
		{
			if(i == len1)
			{
				for(int x = j; x < len2; x++)
				{
					r[k] = r2[x];
					c[k] = c2[x];
					v[k] = v2[x];
					k++;
				}
				break;
			}

			if(j == len2)
			{
				for(int x = i; x < len1; x++)
				{
					r[k] = r1[x];
					c[k] = c1[x];
					v[k] = v1[x];
					k++;
				}
				break;
			}

			if(r1[i] <= r2[j])
			{
				r[k] = r1[i];
				c[k] = c1[i];
				v[k] = v1[i];
				k++;
				i++;
			}
			else
			{
				r[k] = r2[j];
				c[k] = c2[j];
				v[k] = v2[j];
				k++;
				j++;
			}
		}

		delete [] v1; delete [] r1; delete [] c1;
		delete [] v2; delete [] r2; delete [] c2;
	}

	void col_sort()
	/* Sorts elements in each row (after row-sort) by column index. Currently uses insertion sort, which is best for sorting small arrays.
	   Call only after sorting by row! */
	{
		double** vp = new double*[nrows];
		int** cp = new int*[nrows];		// array of pointers to start of each row
		int* lens = new int[nrows];		// length of each row
		for(int i = 0; i < nrows; i++)
			lens[i] = 0;

		// now calculate number of non-zero entries in each row
		int temp=-1, k = -1, iprev = 0;
		for(int i = 0; i < nnz; i++)
		{
			if(temp != row_ind[i])			// when index moves from one row to the next
			{
				k++;
				cp[k] = col_ind + i;
				vp[k] = val + i;
				lens[k]++;
			}
			else lens[k]++;
			temp = row_ind[i];
		}
		/*std::cout << "Lengths of rows: ";
		for(int i = 0; i < nrows; i++)
			std::cout << lens[i] << " ";
		std::cout << std::endl;*/

		// insertion sort
		int tempc = 0; double tempv = 0; int j = 0;
		for(int l = 0; l < nrows; l++)			// for each row
		{
			for(int i = 1; i < lens[l]; i++)
			{
				j = i;
				while(j>0 && cp[l][j] < cp[l][j-1])
				{
					// swap cp[l][j] and cp[l][j-1], same for vp
					tempc = cp[l][j];
					cp[l][j] = cp[l][j-1];
					cp[l][j-1] = tempc;
					tempv = vp[l][j];
					vp[l][j] = vp[l][j-1];
					vp[l][j-1] = tempv;
					j--;
				}
			}
		}

		delete [] vp;
		delete [] cp;
		delete [] lens;
	}

	void sort()
	{
		std::cout << "MatrixCOO: sort(): Sorting the matrix.\n";
		row_sort(val, row_ind, col_ind, nnz);
		//std::cout << "MatrixCOO: sort(): Sorting by column\n";
		col_sort();
		sorted = true;
	}

	void print_data()
	{
		for(int i = 0; i < nnz; i++)
			std::cout << val[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << row_ind[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << col_ind[i] << " ";
		std::cout << std::endl;
	}

	void mprint()
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << get(i,j);
			std::cout << std::endl;
		}
	}

	void fprint(std::ofstream& outfile)
	{
		//outfile << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << get(i,j);
			outfile << std::endl;
		}
	}

	/** Returns product of sparse matrix with vector x and stores it in a.
		NOTE: Parallel version does NOT work!!
	*/
	void multiply(Mat& x, Mat* a, char paralel = 'n')
	{
		//Mat a(x.rows(),x.cols());
		a->zeros();
		//std::cout << "MatrixCOO: multiply(): Rows and columns of a: " << a->rows() << " " << a->cols() << std::endl;
		double temp;
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n"; We don't need our arrays to be sorted.
		#ifdef _OPENMP
		double* val = MatrixCOO::val; int* row_ind = MatrixCOO::row_ind; int* col_ind = MatrixCOO::col_ind; int nnz = MatrixCOO::nnz;
		#endif

		for(int k = 0; k < x.cols(); k++)		// for each column of x
		{
			int i;
			#pragma omp parallel for default(none) private(i,temp) shared(val,row_ind,col_ind,nnz,k,a,x) if(paralel=='y') //num_threads(nthreads_sm)
			for(i = 0; i < nnz; i++)
			{
				temp = val[i]*x.get(col_ind[i],k);
				//(*a)(row_ind[i],k) += val[i]*x(col_ind[i],k);
				//#pragma omp critical (omp_sparse_multiply)
				(*a)(row_ind[i],k) += temp;
			}
		}
		//return a;
	}

	void multiply_parts(Mat* x, Mat* y, Mat* ans, int p)
	/* Like the multiply() method, except the argument matrix is considered x for the first p-1 rows, 0 in the pth row and y for the remaining rows. */
	{
		ans->zeros();
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";
		#ifdef _OPENMP
		double* val = MatrixCOO::val; int* row_ind = MatrixCOO::row_ind; int* col_ind = MatrixCOO::col_ind; int nnz = MatrixCOO::nnz;
		#endif

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i) shared(val,row_ind,col_ind,nnz,k,ans,x,y,p) num_threads(nthreads_sm)
			for(i = 0; i < nnz; i++)
			{
				if(col_ind[i] < p)
					(*ans)(row_ind[i],k) += val[i] * x->get(col_ind[i],k);
				else if(col_ind[i] > p)
					(*ans)(row_ind[i],k) += val[i] * y->get(col_ind[i],k);
			}
		}
	}

	double getelem_multiply_parts(int rownum, Mat* x, Mat* y, int p, double num)
	/* Returns dot product of rownum'th row of this matrix with a certain vector. This vector is composed of x for for the first p-1 elements and y for p+1'th element onwards with num as pth element.
	NOTE: Make sure dimensions of x and y are same. */
	{
		double ans = 0;
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			#ifdef _OPENMP
			double* val = MatrixCOO::val; int* row_ind = MatrixCOO::row_ind; int* col_ind = MatrixCOO::col_ind; int nnz = MatrixCOO::nnz;
			#endif
			int i;
			#pragma omp parallel for default(none) private(i) shared(val,row_ind,col_ind,nnz,k,rownum,p,x,y,num) reduction(+:ans) num_threads(nthreads_sm)
			for(i = 0; i < nnz; i++)
			{
				if(row_ind[i] == rownum)
				{
					if(col_ind[i] < p)
						ans += val[i] * x->get(col_ind[i],k);
					else if(col_ind[i] > p)
						ans += val[i] * y->get(col_ind[i],k);
					else
						ans += val[i] * num;
				}
			}
		}
	}

	void get_diagonal(Mat* D)
	{
		//if(sorted == false) { std::cout << "! MatrixCOO: get_diagonal(): Not sorted!!\n"; return; } does not depend on sorting
		D->zeros();
		for(int i = 0; i < nnz; i++)
		{
			if(row_ind[i] == col_ind[i])
			{
				//(*D)(row_ind[i]) = val[i];
				D->set(row_ind[i],0, val[i]);
			}
		}
	}
	
	void get_lower_triangle(MatrixCOO& L)
	{
		delete [] L.val;
		delete [] L.row_ind;
		delete [] L.col_ind;
		
		L.val = new double[nnz];
		L.row_ind = new int[nnz];
		L.col_ind = new int[nnz];
		L.nrows = nrows;
		L.ncols = ncols;
		
		L.valsize = nnz;
		L.nnz = 0;
		int k = 0;
		
		for(int i = 0; i < nnz; i++)
		{
			if(row_ind[i] > col_ind[i])
			{
				L.val[k] = val[i];
				L.row_ind[k] = row_ind[i];
				L.col_ind[k] = col_ind[i];
				k++;
			}
		}
		
		L.nnz = k;
		
		// shrink to fit
		/*if(L.nnz < nnz)
		{
			for(int i = nnz-1; i >= L.nnz; i--)
			{
				delete L.val+i;
			}
		}*/
	}
	
	void get_upper_triangle(MatrixCOO& L)
	{
		delete [] L.val;
		delete [] L.row_ind;
		delete [] L.col_ind;
		
		L.val = new double[nnz];
		L.row_ind = new int[nnz];
		L.col_ind = new int[nnz];
		L.nrows = nrows;
		L.ncols = ncols;
		
		L.valsize = nnz;
		L.nnz = 0;
		int k = 0;
		
		for(int i = 0; i < nnz; i++)
		{
			if(row_ind[i] < col_ind[i])
			{
				L.val[k] = val[i];
				L.row_ind[k] = row_ind[i];
				L.col_ind[k] = col_ind[i];
				k++;
			}
		}
		
		L.nnz = k;
		
		// shrink to fit
		/*if(L.nnz < nnz)
		{
			for(int i = nnz-1; i >= L.nnz; i--)
			{
				delete L.val+i;
			}
		}*/
	}

	MatrixCOO transpose()
	{
		std::cout << "MatrixCOO: transpose(): Transposing matrix" << std::endl;
		MatrixCOO mat;
		mat.val = new double[nnz];
		mat.row_ind = new int[nnz];
		mat.col_ind = new int[nnz];
		mat.valsize = nnz;
		mat.nrows = nrows; mat.ncols = ncols; mat.sorted = sorted; mat.nnz = nnz;
		int i;

		for(i = 0; i < nnz; i++)
		{
			mat.val[i] = val[i];
			mat.row_ind[i] = col_ind[i];
			mat.col_ind[i] = row_ind[i];
		}

		return mat;
	}

	void combine_sparse_matrices(MatrixCOO& A11, MatrixCOO& A12, MatrixCOO& A21, MatrixCOO& A22)
	/* Combines A, B, C, and D (4 n-by-n matrices) into one 2n-by-2n matrix. CAUTION: All 4 input matrices must have same size! */
	{
		std::cout << "MatrixCOO: combine_sparse_matrices(): Combining matrices" << std::endl;
		//std::cout << "MatrixCOO: combine_sparse_matrices(): A11 dimensions: " << A11.nrows << " " << A11.ncols << std::endl;

		nnz = A11.nnz + A12.nnz + A21.nnz + A22.nnz;
		valsize = nnz+10;
		delete [] val; delete [] row_ind; delete [] col_ind;
		val = new double[valsize]; row_ind = new int[valsize]; col_ind = new int[valsize];
		for(int j = 0; j < valsize; j++)
		{
			row_ind[j] = -1;
			col_ind[j] = -1;
			val[j] = 0.0;
		}
		int i, ii;

		for(i = 0; i < A11.nnz; i++)
		{
			val[i] = A11.val[i];
			row_ind[i] = A11.row_ind[i];
			col_ind[i] = A11.col_ind[i];
		}
		ii = i;

		for(i = ii; i < ii + A12.nnz; i++)
		{
			val[i] = A12.val[i - ii];
			row_ind[i] = A12.row_ind[i - ii];
			col_ind[i] = A12.col_ind[i - ii] + A11.ncols;
		}
		ii = i;

		for(i = ii; i < ii + A21.nnz; i++)
		{
			val[i] = A21.val[i - ii];
			row_ind[i] = A21.row_ind[i-ii] + A11.nrows;
			col_ind[i] = A21.col_ind[i-ii];
		}
		ii = i;

		for(i = ii; i < ii + A22.nnz; i++)
		{
			val[i] = A22.val[i-ii];
			row_ind[i] = A22.row_ind[i-ii] + A11.nrows;
			col_ind[i] = A22.col_ind[i-ii] + A11.ncols;
		}
		ii = i;

		for(i = ii; i < valsize; i++)
		{
			val[i] = 0;
			row_ind[i] = -1;
			col_ind[i] = -1;
		}
	}
};

class MatrixCRS : public SparseMatrix
{
	int* col_ind;
	int* row_ptr;
	int nnz;
	int nrows;
	int ncols;
	int valsize;		// total size allocated to array val, and therefore also to array col_ind

public:
	MatrixCRS(int num_rows, int num_cols) : SparseMatrix(num_rows, num_cols)
	{
		val = new double[num_rows*2];
		col_ind = new int[num_rows*2];
		valsize = num_rows*2;
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0.0;
			col_ind[i] = -1;
		}
		row_ptr = new int[nrows+1];
		for(int i = 0; i < nrows; i++) row_ptr[i] = -1;
		row_ptr[nrows+1] = 0;
		nnz = 0;
	}
	~MatrixCRS()
	{
		delete [] val;
		delete [] col_ind;
		delete [] row_ptr;
	}

	void set(int x, int y, double value)
	{
		if(col_ind[valsize-1] != -1)		// re-allocate if no space is available
		{
			double* oldval = val;
			int* oldcol = col_ind;
			val = new double[2*valsize];
			col_ind = new int[2*valsize];
			for(int i = 0; i < valsize; i++)		// copy data over
			{
				val[i] = oldval[i];
				col_ind[i] = oldcol[i];
			}
			for(int i = valsize; i < 2*valsize; i++)		// initialize extra space with dummy data
			{
				val[i] = 0;
				col_ind[i] = -1;
			}
			delete [] oldval;
			delete [] oldcol;
		}

		// if value is not zero

		int ind = row_ptr[x];
		int nind = row_ptr[x+1];
		if(ind == nind)				// only zeros in this row
		{
			for(int i = nind; i < nnz; i++)
			{
				val[i+1] = val[i];				// create space for this element
				col_ind[i+1] = col_ind[i];
			}
		}
	}

	double get(int x, int y=0)
	{
		int ind = row_ptr[x];
		//if(col_ind[ind] > y) return 0.0;
		int nind = row_ptr[x+1];
		for(int i = ind; i <= nind-1; i++)
		{
			if(col_ind[i]==y) return val[i];
		}
		return 0.0;
	}

	double& operator()(int x, int y=0)
	{

	}
};

}
