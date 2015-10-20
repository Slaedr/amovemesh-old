/* Classes to reconstruct a C^1 (or C^2) piecewise polynomial (cubic) boundary from a piecewise linear boundary given by a linear mesh.
   Aditya Kashi
   September 4, 2015
*/

#ifndef __AMESH2DGENERAL_H
#include "amesh2d.hpp"
#endif

using namespace std;
using namespace amat;
using namespace acfd;

namespace acfd {

typedef double(*Basispointer)(double);

/* Class HermiteSpline2d constructs one spline curve from all the boundary faces of a mesh which have markers contained in rfl.
   Currently reconstructs a C^1 boundary.
   NOTE: Make sure the boundaries with markers in rfl are contiguous.
*/
class HermiteSpline2d
{
	typedef double(HermiteSpline2d::*Basispointer)(double);
	UMesh2d* m;
	Matrix<int> rfl;										// the markers of the boundaries to be reconstructed
	int nseg;												// number of curve segments - equal to number of boundary faces to reconstruct
	int ndeg;												// degree of spline interpolation - usually 3
	Matrix<double> cf;										// to store geometric coefficients
	Matrix<int> bfaces;										// boundary face data copied from m->bface
	Matrix<int> obfaces;									// ordered bfaces
	Matrix<double> gallfa;									// Normals etc of boundary faces in bfaces
	Matrix<int> facep;										// Maps mesh bface number to the corresponding data in bfaces and gallfa -
															//  will help in getting spline data based on mesh's original bface numbering.
	bool isclosed;											// True if the faces to be reconstructed for a closed boundary.
	Basispointer* F;										// Array of function pointers for Hermite basis functions

	// definitions of the 4 Hermite basis functions
	double f0(double u) { return 2*pow(u,3) - 3*pow(u,2)+1; }
	double f1(double u) { return -2*pow(u,3) + 3*pow(u,2); }
	double f2(double u) { return pow(u,3) - 2*u*u + u; }
	double f3(double u) { return u*u*u - u*u; }

public:
	void setup(UMesh2d* mesh, Matrix<int> rflags)
	{
		m = mesh;
		rfl = rflags;
		ndeg = 3;

		F = new Basispointer[ndeg+1];

		F[0] = &HermiteSpline2d::f0;
		F[1] = &HermiteSpline2d::f1;
		F[2] = &HermiteSpline2d::f2;
		F[3] = &HermiteSpline2d::f3;
		//cout << (this->*F[0])(2.1);

		facep.setup(m->gnface(),1);

		nseg = 0;
		for(int ifa = 0; ifa < m->gnface(); ifa++)
		{
			for(int imark = 0; imark < rfl.msize(); imark++)
				if(m->gbface(ifa,m->gnnofa()) == rfl(imark)) nseg++;
		}

		cf.setup(nseg,ndeg+1);
		bfaces.setup(nseg,m->gnnofa());

		nseg = 0;
		for(int ifa = 0; ifa < m->gnface(); ifa++)
		{
			for(int imark = 0; imark < rfl.msize(); imark++)
				if(m->gbface(ifa,m->gnnofa()) == rfl(imark))
				{
					for(int j = 0; j < m->gnnofa(); j++)
						bfaces(nseg,j) = m->gbface(ifa,j);
					facep(ifa) = nseg;
					nseg++;
					break;
				}
		}

		isclosed = true;
		// now re-order bfaces such that the boundary faces are contiguous, and point 2 of face i is point 1 of face i+1
		// this is so that it's easy to get consistent tangents
	}

	~HermiteSpline2d()
	{
		delete [] F;
	}

	void compute_splines()
	{
		// iterate over boundary faces of the mesh
		for(int iseg = 0; iseg < nseg; iseg++)
		{
			// get tangent for this face
			// if edge is bounded by (x1,y1) and (x2,y2), tangent to edge is (x2-x1)i + (y2-y1)j
			// But we need to ensure that out of the two poosible tangents, the correct one is chosen for consistency.
		}
	}
};

}
