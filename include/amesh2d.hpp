// Data structure and setup for 2D unstructured mesh.
// Aditya Kashi
// Aug 12, 2015

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif
#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif
#ifndef _GLIBCXX_STRING
#include <string>
#endif
#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif
#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif
#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif
#ifndef __ADATASTRUCTURES_H
#include <adatastructures.hpp>
#endif

#define __AMESH2DGENERAL_H

using namespace std;
using namespace amat;

namespace acfd {

const double PI = 3.14159265358979323846;

/** Class UMesh2d is a general mesh class for unstructured mesh (with 1 kind of element throughout - hybrid mesh is NOT supporte) */
class UMesh2d
{
private:
	int npoin;		//< Number of nodes
	int nelem;		//< Number of elements
	int nface;		//< Number of boundary faces
	int ndim;		//< Dimension of the mesh
	int nnode;		//< number of nodes to an element
	int nfael;		//< number of faces to an element (equal to number of edges to an element in 2D)
	int nnofa;		//< number of node in a face -- needs to be generalized in case of general grids
	int naface;		//< total number of (internal and boundary) faces
	int nbface;		//< number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	int nbpoin;		//< number of boundary points
	int nbtag;		//< number of tags for each boundary face
	int ndtag;		//< number of tags for each element
	Matrix<double> coords;			//< Specifies coordinates of each node
	Matrix<int> inpoel;				//< Interconnectivity matrix: lists node numbers of nodes in each element
	Matrix<int> bface;				//< Boundary face data: lists nodes belonging to a boundary face and contains boudnary markers
	Matrix<double> vol_regions;		//< to hold volume region markers, if any

	Matrix<int> esup_p;
	Matrix<int> esup;
	Matrix<int> psup_p;
	Matrix<int> psup;
	Matrix<int> esuel;
	Matrix<int> intfac;
	Matrix<int> intfacbtags;		//< to hold boundary tags corresponding to intfac

	Matrix<int> bpoints;
	//< bpoints contains: bpoints(0) = global point number, bpoints(1) = first containing intfac face (face with intfac's second point as this point),
	//<  bpoints(2) = second containing intfac face (face with intfac's first point as this point)

	Matrix<int> bpointsb;
	//< Like bpoints, but stores bface numbers corresponding to each face, rather than intfac faces

	Matrix<int> bfacebp;
	//< Stores boundary-points numbers (defined by bpointsb) of the two points making up a particular bface.

	Matrix<int> bifmap;				//< relates boundary faces in intfac with bface, ie, bifmap(intfac no.) = bface no.
	Matrix<int> ifbmap;				//< relates boundary faces in bface with intfac, ie, ifbmap(bface no.) = intfac no.
	bool isBoundaryMaps;			//< Specifies whether bface-intfac maps have been created

	bool alloc_jacobians;
	Matrix<double> jacobians;		//< Contains jacobians of each (linear) element

public:

	UMesh2d() { alloc_jacobians = false; }

	UMesh2d(const UMesh2d& other)
	{
		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		nbtag = other.nbtag;
		ndtag = other.ndtag;
		nbpoin = other.nbpoin;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		vol_regions = other.vol_regions;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		intfacbtags = other.intfacbtags;
		bpoints = other.bpoints;
		bpointsb = other.bpointsb;
		bfacebp = other.bfacebp;
		bifmap = other.bifmap;
		ifbmap = other.ifbmap;
		isBoundaryMaps = other.isBoundaryMaps;
		//gallfa = other.gallfa;
		alloc_jacobians = other.alloc_jacobians;
		jacobians = other.jacobians;
	}

	UMesh2d& operator=(const UMesh2d& other)
	{
		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		nbtag = other.nbtag;
		ndtag = other.ndtag;
		nbpoin = other.nbpoin;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		vol_regions = other.vol_regions;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		intfacbtags = other.intfacbtags;
		bpoints = other.bpoints;
		bpointsb = other.bpointsb;
		bfacebp = other.bfacebp;
		bifmap = other.bifmap;
		ifbmap = other.ifbmap;
		isBoundaryMaps = other.isBoundaryMaps;
		//gallfa = other.gallfa;
		alloc_jacobians = other.alloc_jacobians;
		jacobians = other.jacobians;
		return *this;
	}

	double gcoords(int pointno, int dim) const
	{
		return coords.get(pointno,dim);
	}
	int ginpoel(int elemno, int locnode) const
	{
		return inpoel.get(elemno, locnode);
	}
	int gbface(int faceno, int val) const
	{
		return bface.get(faceno, val);
	}

	void setcoords(Matrix<double>* c)
	{ coords = *c; }

	void setinpoel(Matrix<int>* inp)
	{ inpoel = *inp; }

	void setbface(Matrix<int>* bf)
	{ bface = *bf; }

	void modify_bface_marker(int iface, int pos, int number)
	{ bface(iface, pos) = number; }

	Matrix<double>* getcoords()
	{ return &coords; }

	int gesup(int i) { return esup.get(i); }
	int gesup_p(int i) { return esup_p.get(i); }
	int gpsup(int i) { return psup.get(i); }
	int gpsup_p(int i) { return psup_p.get(i); }
	int gesuel(int ielem, int jnode) { return esuel.get(ielem, jnode); }
	int gintfac(int face, int i) const { return intfac.get(face,i); }
	int gintfacbtags(int face, int i) { return intfacbtags.get(face,i); }
	int gbpoints(int poin, int i) { return bpoints.get(poin,i); }
	int gbpointsb(int poin, int i) { return bpointsb.get(poin,i); }
	int gbfacebp(int iface, int i) { return bfacebp.get(iface,i); }
	int gbifmap(int intfacno) { return bifmap.get(intfacno); }
	int gifbmap(int bfaceno) { return ifbmap.get(bfaceno); }
	double gjacobians(int ielem) { return jacobians(ielem,0); }

	int gnpoin() { return npoin; }
	int gnelem() { return nelem; }
	int gnface() { return nface; }
	int gnbface() { return nbface; }
	int gnnode() { return nnode; }
	int gndim() { return ndim; }
	int gnaface() {return naface; }
	int gnfael() { return nfael; }
	int gnnofa() { return nnofa; }
	int gnbtag() {return nbtag; }
	int gndtag() { return ndtag; }
	int gnbpoin() { return nbpoin; }

	void readDomn(string mfile)
	{
		/** Reads Professor Luo's mesh file, which I call the 'domn' format.
		   NOTE: Make sure nfael and nnofa are mentioned after ndim and nnode in the mesh file.
		*/

		ifstream infile(mfile);
		
		// Do file handling here to populate npoin and nelem
		cout << "UTriMesh: Reading mesh file...\n";
		char ch = '\0'; int dum = 0; double dummy;

		infile >> dum;
		infile >> ch;
		for(int i = 0; i < 4; i++)		//skip 4 lines
			do
				ch = infile.get();
			while(ch != '\n');
		infile >> ndim;
		infile >> nnode;
		infile >> nfael;
		infile >> nnofa;
		infile >> ch;			//get the newline
		do
			ch = infile.get();
		while(ch != '\n');
		infile >> nelem; infile >> npoin; infile >> nface;
		infile >> dummy; 				// get time
		ch = infile.get();			// clear newline

		cout << "UTriMesh: Number of elements: " << nelem << ", number of points: " << npoin << ", number of nodes per element: " << nnode << endl;
		cout << "Number of boundary faces: " << nface << ", Number of dimensions: " << ndim;

		nbtag = 2;
		ndtag = 2;

		//cout << "\nUTriMesh: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		//cout << "\nUTriMesh: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		//cout << "UTriMesh: Allocating bface...\n";
		bface.setup(nface, nnofa + nbtag, ROWMAJOR);
		
		//cout << "UTriMesh: Allocation done.";

		do
			ch = infile.get();
		while(ch != '\n');

		//now populate inpoel
		for(int i = 0; i < nelem; i++)
		{
			infile >> dum;

			for(int j = 0; j < nnode; j++)
				infile >> inpoel(i,j);

			do
				ch = infile.get();
			while(ch != '\n');
		}
		cout << "\nUTriMesh: Populated inpoel.";

		//Correct inpoel:
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
				inpoel(i,j)--;
		}

		ch = infile.get();
		do
			ch = infile.get();
		while(ch != '\n');

		// populate coords
		for(int i = 0; i < npoin; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim; j++)
				infile >> coords(i,j);
		}
		cout << "\nUTriMesh: Populated coords.\n";

		// skip initial conditions
		ch = infile.get();
		for(int i = 0; i < npoin+2; i++)
		{
			do
				ch = infile.get();
			while(ch != '\n');
		}
		
		// populate bface
		for(int i = 0; i < nface; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim + nbtag; j++)
			{
				infile >> bface(i,j);
			}
			if (i==nface-1) break;
			do
				ch = infile.get();
			while(ch!='\n');
		}
		cout << "UTriMesh: Populated bface. Done reading mesh.\n";
		//correct first 2 columns of bface
		for(int i = 0; i < nface; i++)
			for(int j = 0; j < 2; j++)
				bface(i,j)--;

		infile.close();

		vol_regions.setup(nelem, ndtag);
		vol_regions.zeros();
	}

	/// Reads mesh from Gmsh 2 format file
	void readGmsh2(string mfile, int dimensions)
	{
		cout << "UMesh2d: readGmsh2(): Reading mesh file...\n";
		int dum; double dummy; string dums; char ch;
		ndim = dimensions;

		ifstream infile(mfile);
		for(int i = 0; i < 4; i++)		//skip 4 lines
			do
				ch = infile.get();
			while(ch != '\n');

		infile >> npoin;
		cout << "UMesh2d: readGmsh2(): No. of points = " << npoin << endl;
		coords.setup(npoin,ndim);

		// read coords of points
		for(int i = 0; i < npoin; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim; j++)
				infile >> coords(i,j);
			if(ndim < 3) infile >> dummy;
		}
		infile >> dums;
		infile >> dums;
		//cout << "UMesh2d: readGmsh2(): coords read." << endl;

		int nelm, elmtype, nbtags, ntags;
		/// elmtype is the standard element type in the Gmsh 2 mesh format - of either faces or elements
		ndtag = 0; nbtag = 0;
		infile >> nelm;
		Matrix<int> elms(nelm,20);
		nface = 0; nelem = 0;
		//cout << "UMesh2d: readGmsh2(): Total number of elms is " << nelm << endl;

		for(int i = 0; i < nelm; i++)
		{
			infile >> dum;
			infile >> elmtype;
			/// assumption here is that elmtype is same for all faces and same for all elements
			switch(elmtype)
			{
				case(1): // linear edge
					nnofa = 2;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(8): // quadratic edge
					nnofa = 3;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(2): // linear triangles
					nnode = 3;
					nfael = 3;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(3):	// linear quads
					nnode = 4;
					nfael = 4;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(9):	// quadratic triangles
					nnode = 6;
					nfael = 3;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(16):	// quadratic quad (8 nodes)
					nnode = 8;
					nfael = 4;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(10):	// quadratic quad (9 nodes)
					nnode = 9;
					nfael = 4;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				default:
					cout << "! UMesh2d: readGmsh2(): Element type not recognized. Setting as linear triangle." << endl;
					nnode = 3;
					nfael = 3;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
			}
		}
		//cout << "UMesh2d: readGmsh2(): Done reading elms" << endl;

		if(nface > 0)
			bface.setup(nface, nnofa+nbtag);
		else cout << "UMesh2d: readGmsh2(): NOTE: There is no boundary data!" << endl;

		inpoel.setup(nelem, nnode);
		vol_regions.setup(nelem, ndtag);

		cout << "UMesh2d: readGmsh2(): Done. No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ",\n number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << endl;

		// write into inpoel and bface
		// the first nface rows to be read are boundary faces
		for(int i = 0; i < nface; i++)
		{
			for(int j = 0; j < nnofa; j++)
				bface(i,j) = elms(i,j)-1;			// -1 to correct for the fact that our numbering starts from zero
			for(int j = nnofa; j < nnofa+nbtag; j++)
				bface(i,j) = elms(i,j);
		}
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
				inpoel(i,j) = elms(i+nface,j)-1;
			for(int j = 0; j < ndtag; j++)
				vol_regions(i,j) = elms(i+nface,j+nnode);
		}
		infile.close();
	}

	void readYibinMesh(string meshn)
	{
		/** Reads Yibin Wang's dat mesh format, which does not contain any boundary information.
			Assumes various properties of the mesh, as the file does not specify them explicitly.
			Also carries out a process to determine boundary data bface and bpoints, which only works
			for structured meshes.
		*/

		ifstream fin(meshn);

		fin >> npoin >> nelem;

		// some assumptions for properties not specified in the file
		ndim = 2; nnode = 4; nnofa = 2; nfael = 4; nbtag = 2; ndtag = 2;
		int a;

		coords.setup(npoin,ndim);
		inpoel.setup(nelem,nnode);
		vol_regions.setup(nelem, ndtag);

		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			for(int idim = 0; idim < ndim; idim++)
				fin >> coords(ipoin,idim);
		}
		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int inode = 0; inode < nnode; inode++) {
				fin >> a;
				inpoel(ielem,inode) = a - 1;
			}
		}

		// finished reading mesh file
		fin.close();

		// Now populate bface and bpoints
		// NOTE: only works for structured mesh!

		Matrix<int> lpoin(npoin,1);
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int inode = 0; inode < nnode; inode++)
				lpoin(inpoel(ielem,inode))++;
		}

		nbpoin = 0;
		nface = 0;
		for(int i = 0; i < npoin; i++)
			if(lpoin(i) <= 3) nbpoin++;

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			int nextnode;
			for(int inode = 0; inode < nnode; inode++)
			{
				nextnode = perm(0,nnode-1,inode,1);
				if(lpoin.get(inpoel(ielem,inode)) <= 3 && lpoin.get(inpoel.get(ielem,nextnode)) <= 3) {
					nface++;
				}
			}
		}

		bface.setup(nface,nnofa+nbtag);
		bface.zeros();
		Matrix<int> trav(npoin,1);
		trav.zeros();
		int k = 0, l = 0;

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			int thisgnode, nextnode, nextgnode;
			for(int inode = 0; inode < nnode; inode++)
			{
				thisgnode = inpoel(ielem,inode);
				nextnode = perm(0,nnode-1,inode,1);
				nextgnode = inpoel(ielem,nextnode);
				if(lpoin.get(thisgnode) <= 3 && lpoin.get(nextgnode) <= 3) {
					bface(k,0) = nextgnode;
					bface(k,1) = thisgnode;
					bface(k,2) = 1;
					bface(k,3) = 0;
					k++;
				}
			}
		}

		// we now have bpoints and bface for a structured mesh.
		cout << "UTriMesh: readYibinMesh(): " << npoin << " points, " << nelem << " elements, " << nface << " boundary faces, " << nbpoin << " boundary points." << endl;

		compute_boundary_points();
		
		// assign boundary markers
		int face1 = 0, face2 = 125, curface, nextface, nextpoin;
		nextface = face1;
		while(true)
		{
			curface = nextface;
			bface(curface,nnofa) = 1;
			nextpoin = bfacebp(curface,1);
			nextface = bpointsb(nextpoin,2);
			if(nextface == face1) break;
		}

		nextface = face2;
		while(true)
		{
			curface = nextface;
			bface(curface,nnofa) = 2;
			nextpoin = bfacebp(curface,1);
			nextface = bpointsb(nextpoin,2);
			if(nextface == face2) break;
		}
	}

	void compute_boundary_points()
	/** Stores for each boundary point: the associated global point number and the two bfaces associated with it.
	*/
	{
		cout << "UMesh2d: compute_boundary_points(): Calculating bpointsb structure"<< endl;

		// first, get number of boundary points

		nbpoin = 0;
		Matrix<int> flagb(npoin,1);
		flagb.zeros();
		for(int iface = 0; iface < nface; iface++)
		{
			for(int inofa = 0; inofa < nnofa; inofa++)
				flagb(bface(iface,inofa)) = 1;
		}
		for(int ipoin = 0; ipoin < npoin; ipoin++)
			nbpoin += flagb(ipoin);

		cout << "UMesh2d: compute_boundary_points(): No. of boundary points = " << nbpoin << endl;

		bpointsb.setup(nbpoin,3);
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < 3; j++)
				bpointsb(i,j) = -1;

		bfacebp.setup(nface,nnofa);
		
		Matrix<double> lpoin(npoin,1);

		int bp = 0;

		// Next, populate bpointsb by iterating over faces. Also populate bfacebp, which holds the boundary points numbers of the 2 points in a bface.
		lpoin.zeros();		// lpoin will be 1 if the point has been visited
		for(int iface = 0; iface < nface; iface++)
		{
			int p1, p2;
			p1 = bface(iface,0);
			p2 = bface(iface,1);

			if(lpoin(p1) == 0)	// if this point has not been visited before
			{
				bpointsb(bp,0) = p1;
				bpointsb(bp,2) = iface;
				bfacebp(iface,0) = bp;
				bp++;
				lpoin(p1) = 1;
			}
			else
			{
				// search bpoints for point p1
				int ibp=-1;
				for(int i = 0; i < bp; i++)
				{
					if(bpointsb(i,0) == p1) ibp = i;
				}
				if(ibp==-1) cout << "UMesh2d: compute_boundary_points(): Point not found!" << endl;
				bpointsb(ibp,2) = iface;
				bfacebp(iface,0) = ibp;
			}

			if(lpoin(p2) == 0)	// if this point has not been visited before
			{
				bpointsb(bp,0) = p2;
				bpointsb(bp,1) = iface;
				bfacebp(iface,1) = bp;
				bp++;
				lpoin(p2) = 1;
			}
			else
			{
				// search bpoints for point p2
				int ibp=-1;
				for(int i = 0; i < bp; i++)
				{
					if(bpointsb(i,0) == p2) ibp = i;
				}
				if(ibp==-1) cout << "UMesh2d: compute_boundary_points(): Point not found!" << endl;
				bpointsb(ibp,1) = iface;
				bfacebp(iface,1) = ibp;
			}
		}
	}

	void printmeshstats()
	{
		cout << "UMesh2d: No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ", number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << endl;
	}

	void writeGmsh2(string mfile)
	{
		cout << "UMesh2d: writeGmsh2(): writing mesh to file " << mfile << endl;
		// decide element type first, based on nfael/nnode and nnofa
		int elm_type = 2;
		int face_type = 1;
		if(nnode == 4)
			elm_type = 3;
		else if(nnode == 6)
			elm_type = 9;
		else if(nnode == 8)
			elm_type = 16;
		else if(nnode==9)
			elm_type = 10;

		if(nnofa == 3)
			face_type = 8;

		ofstream outf(mfile);
		//cout << "nodes\n";
		outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		outf << "$Nodes\n" << npoin << '\n';
		for(int ip = 0; ip < npoin; ip++)
		{
			outf << ip+1;
			for(int j = 0; j < ndim; j++)
			 	outf << " " << coords(ip,j);
			for(int j = 3-ndim; j > 0; j--)
				outf << " " << 0.0;
			outf << '\n';
		}
		outf << "$EndNodes\n";

		//cout << "boundary faces\n";
		outf << "$Elements\n" << nelem+nface << '\n';
		// boundary faces first
		for(int iface = 0; iface < nface; iface++)
		{
			outf << iface+1 << " " << face_type << " " << nbtag;
			for(int i = nnofa; i < nnofa+nbtag; i++)
				outf << " " << bface(iface,i);			// write tags
			for(int i = 0; i < nnofa; i++)
				outf << " " << bface(iface,i)+1;		// write nodes
			outf << '\n';
		}
		//cout << "elements\n";
		for(int iel = 0; iel < nelem; iel++)
		{
			outf << nface+iel+1 << " " << elm_type << " " << ndtag;
			for(int i = 0; i < ndtag; i++)
				outf << " " << vol_regions(iel,i);
			for(int i = 0; i < nnode; i++)
				outf << " " << inpoel(iel,i)+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
	}

	void compute_jacobians()
	{
		if(nnode == 3 || nnode == 4)
		{
			if (alloc_jacobians == false)
			{
				jacobians.setup(nelem, 1);
				alloc_jacobians = true;
			}

			for(int i = 0; i < gnelem(); i++)
			{
				jacobians(i,0) = gcoords(ginpoel(i,0),0)*(gcoords(ginpoel(i,1),1) - gcoords(ginpoel(i,2),1)) - gcoords(ginpoel(i,0),1)*(gcoords(ginpoel(i,1),0)-gcoords(ginpoel(i,2),0)) + gcoords(ginpoel(i,1),0)*gcoords(ginpoel(i,2),1) - gcoords(ginpoel(i,2),0)*gcoords(ginpoel(i,1),1);
			}
		}
		else {
			cout << "UMesh2d: compute_jacobians(): ! Mesh is not linear. Cannot compute jacobians." << endl;
		}
	}

	void detect_negative_jacobians(ofstream& out)
	{
		bool flagj = false;
		int nneg = 0;
		for(int i = 0; i < nelem; i++)
		{
			if(jacobians(i,0) <= 1e-15) {
				out << i << " " << jacobians(i,0) << '\n';
				flagj = true;
				nneg++;
			}
		}
		if(flagj == true) cout << "UMesh2d: detect_negative_jacobians(): There exist " << nneg << " element(s) with negative jacobian!!\n";
	}

	void compute_topological()
	{
		/// Computes data structures for elements surrounding point (esup), points surrounding point (psup), elements surrounding elements (esuel),
		/// elements surrounding faces along with points in faces (intfac), and also
		/// a list of boundary points with correspong global point numbers and containing boundary faces (according to intfac) (bpoints).

		/// NOTE: (1) Use only after setup()
		///		  (2) Currently only works for linear mesh

		cout << "UMesh2d: compute_topological(): Calculating and storing topological information...\n";
		/// 1. Elements surrounding points
		//cout << "UMesh2d: compute_topological(): Elements surrounding points\n";
		esup_p.setup(npoin+1,1,ROWMAJOR);
		esup_p.zeros();

		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				esup_p(inpoel(i,j)+1,0) += 1;	// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
			}
		}
		// Now make the members of esup_p cumulative
		for(int i = 1; i < npoin+1; i++)
			esup_p(i,0) += esup_p(i-1,0);
		// Now populate esup
		esup.setup(esup_p(npoin,0),1,ROWMAJOR);
		esup.zeros();
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				int ipoin = inpoel(i,j);
				esup(esup_p(ipoin,0),0) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
				esup_p(ipoin,0) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
			}
		}
		//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
		// So correct this.
		for(int i = npoin; i >= 1; i--)
			esup_p(i,0) = esup_p(i-1,0);
		esup_p(0,0) = 0;
		// Elements surrounding points is now done.

		/// 2. Points surrounding points
		cout << "UMesh2d: compute_topological(): Points surrounding points\n";
		psup_p.setup(npoin+1,1,ROWMAJOR);
		psup_p.zeros();
		psup_p(0,0) = 0;
		Matrix<int> lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize this vector to -1
		int istor = 0;

		// first pass: calculate storage needed for psup
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number

				// find local node number of ip in ielem
				int inode;
				for(int jnode = 0; jnode < nnode; jnode++)
					if(inpoel(ielem,jnode) == ip) inode = jnode;

				vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
				for(int j = 0; j < nnode; j++)
					nbd[j] = false;

				if(nnode == 3)
					for(int i = 0; i < nbd.size(); i++)
						nbd[i] = true;
				else if(nnode == 4)
					for(int jnode = 0; jnode < nnode; jnode++)
					{
						if(jnode == perm(0,nnode-1,inode,1) || jnode == perm(0,nnode-1, inode, -1))
							nbd[jnode] = true;
					}

				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
					{
						istor++;
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
			psup_p(ip+1,0) = istor;
		}

		psup.setup(istor,1,ROWMAJOR);
		//cout << "+++ " << istor << endl;

		//second pass: populate psup
		istor = 0;
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize lpoin to -1
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number

				// find local node number of ip in ielem
				int inode;
				for(int jnode = 0; jnode < nnode; jnode++)
					if(inpoel(ielem,jnode) == ip) inode = jnode;

				vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
				for(int j = 0; j < nnode; j++)
					nbd[j] = false;

				if(nnode == 3)
					for(int i = 0; i < nbd.size(); i++)
						nbd[i] = true;
				else if(nnode == 4)
					for(int jnode = 0; jnode < nnode; jnode++)
					{
						if(jnode == perm(0,nnode-1,inode,1) || jnode == perm(0,nnode-1, inode, -1))
							nbd[jnode] = true;
					}

				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
					{
						psup(istor,0) = jpoin;
						istor++;
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
		}
		//Points surrounding points is now done.

		/// 3. Elements surrounding elements
		//cout << "UMesh2d: compute_topological(): Elements surrounding elements...\n";

		//Matrix<int> lpoin(npoin,1);
		esuel.setup(nelem, nfael, ROWMAJOR);
		for(int ii = 0; ii < nelem; ii++)
			for(int jj = 0; jj < nfael; jj++)
				esuel(ii,jj) = -1;
		Matrix<int> lpofa(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])
		/*lpofa(0,0) = 1; lpofa(0,1) = 2;
		lpofa(1,0) = 2; lpofa(1,1) = 0;
		lpofa(2,0) = 0; lpofa(2,1) = 1;*/
		for(int i = 0; i < nfael; i++)
		{
			for(int j = 0; j < nnofa; j++)
			{
				lpofa(i,j) = perm(0,nnode-1,i,j);
			}
		}
		//lpofa.mprint();
		Matrix<int> lhelp(nnofa,1);
		lhelp.zeros();
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int ifael = 0; ifael < nfael; ifael++)
			{
				for(int i = 0; i < nnofa; i++)
				{
					lhelp(i,0) = inpoel(ielem, lpofa(ifael,i));	// lhelp stores global node nos. of current face of current element
					lpoin(lhelp(i,0)) = 1;
				}
				int ipoin = lhelp(0);
				for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
				{
					int jelem = esup(istor);
					if(jelem != ielem)
					{
						for(int jfael = 0; jfael < nfael; jfael++)
						{
							//Assume that no. of nodes in face ifael is same as that in face jfael
							int icoun = 0;
							for(int jnofa = 0; jnofa < nnofa; jnofa++)
							{
								int jpoin = inpoel(jelem, lpofa(jfael,jnofa));
								if(lpoin(jpoin)==1) icoun++;
							}
							if(icoun == nnofa)		// nnofa is 2
							{
								esuel(ielem,ifael) = jelem;
								esuel(jelem,jfael) = ielem;
							}
						}
					}
				}
				for(int i = 0; i < nnofa; i++)
					lpoin(lhelp(i)) = 0;
			}
		}

		/** Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac. Also computes unit normals to, and lengths of, each face as well as boundary flags of boundary faces, in gallfa.
		The orientation of the face is such that the element with smaller index is always to the left of the face, while the element with greater index is always to the right of the face.
		NOTE: After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.*/

		//cout << "UMesh2d: compute_topological(): Computing intfac..." << endl;
		nbface = naface = 0;
		// first run: calculate nbface
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				//int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je == -1)
				{
					//esuel(ie,in) = nelem+nbface;
					nbface++;
				}
			}
		}
		cout << "UMesh2d: compute_topological(): Number of boundary faces = " << nbface << endl;
		// calculate number of internal faces
		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				//int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je > ie && je < nelem) naface++;
			}
		}
		cout << "UMesh2d: compute_topological(): Number of all faces = " << naface << endl;

		//allocate intfac
		intfac.setup(naface,nnofa+2,ROWMAJOR);

		//reset face totals
		nbface = naface = 0;

		//second run: populate intfac
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je == -1)
				{
					esuel(ie,in) = nelem+nbface;
					intfac(nbface,0) = ie;
					intfac(nbface,1) = nelem+nbface;
					intfac(nbface,2) = inpoel(ie,in);
					intfac(nbface,3) = inpoel(ie,in1);

					nbface++;
				}
			}
		}
		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je > ie && je < nelem)
				{
					intfac(naface,0) = ie;
					intfac(naface,1) = je;
					intfac(naface,2) = inpoel(ie,in);
					intfac(naface,3) = inpoel(ie,in1);
					naface++;
				}
			}
		}

		/// Finally, calculates bpoints.

		//first get number of bpoints
		nbpoin = 0;
		Matrix<int> isbpflag(npoin,1);
		isbpflag.zeros();
		for(int i = 0; i < nface; i++)
		{
			for(int j = 0; j < nnofa; j++)
				isbpflag(bface(i,j)) = 1;
		}
		for(int i = 0; i < npoin; i++)
			if(isbpflag(i)==1) nbpoin++;

		cout << "UMesh2d: compute_topological(): Number of boundary points = " << nbpoin << endl;

		//Allocate bpoints
		bpoints.setup(nbpoin,3);		// We need 1 field for global point number and in 2D linear meshes, we need 2 more for surrounding faces
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < 3; j++)
				bpoints(i,j) = -1;

		int bp = 0;

		// Next, populate bpoints by iterating over intfac faces
		lpoin.zeros();		// lpoin will be 1 if the point has been visited
		for(int iface = 0; iface < nbface; iface++)
		{
			int p1, p2;
			p1 = intfac(iface,2+0);
			p2 = intfac(iface,2+1);

			if(lpoin(p1) == 0)	// if this point has not been visited before
			{
				bpoints(bp,0) = p1;
				bpoints(bp,2) = iface;
				bp++;
				lpoin(p1) = 1;
			}
			else
			{
				// search bpoints for point p1
				int ibp=-1;
				for(int i = 0; i < bp; i++)
				{
					if(bpoints(i,0) == p1) ibp = i;
				}

				bpoints(ibp,2) = iface;
			}

			if(lpoin(p2) == 0)	// if this point has not been visited before
			{
				bpoints(bp,0) = p2;
				bpoints(bp,1) = iface;
				bp++;
				lpoin(p2) = 1;
			}
			else
			{
				// search bpoints for point p2
				int ibp;
				for(int i = 0; i < bp; i++)
				{
					if(bpoints(i,0) == p2) ibp = i;
				}

				bpoints(ibp,1) = iface;
			}
		}
	}

	void compute_boundary_maps()
	{
		// iterate over bfaces and find corresponding intfac face for each bface
		bifmap.setup(nbface,1);
		ifbmap.setup(nbface,1);

		vector<int> fpo(nnofa);

		for(int ibface = 0; ibface < nface; ibface++)
		{
			for(int i = 0; i < nnofa; i++)
				fpo[i] = bface(ibface,i);

			int inface = -1;

			// iterate over intfacs - check if bface ibface matches intfac iface, for each iface
			for(int iface = 0; iface < nbface; iface++)
			{
				bool final1 = true;

				vector<bool> inter(nnofa);
				for(int b = 0; b < nnofa; b++)
					inter[b] = false;						// initially set all bools to false

				for(int j = 0; j < nnofa; j++)
				{
					for(int k = 0; k < nnofa; k++)
						if(fpo[j] == intfac(iface, 2+k)) {
							inter[j] = true;			// if jth node of ibface has a node of iface, it belongs to iface; set the corresp. boolean to true
							break;
						}
				}

				/*for(int i = 0; i < nnofa; i++)
					cout << inter[i];
				cout << endl;*/

				for(int b = 0; b < nnofa; b++)
					if(inter[b] == false) final1 = false;						// if any node of ibface failed to find a node of iface, ibface is not the same as iface

				if(final1 == true) inface = iface;
			}

			if(inface != -1) {
				bifmap(inface) = ibface;
				ifbmap(ibface) = inface;
			}
			else {
				cout << "UMesh2d: compute_boundary_maps(): ! intfac face corresponding to " << ibface << "th bface not found!!" << endl;
				continue;
			}
		}
		isBoundaryMaps = true;
	}

	void writeBoundaryMapsToFile(string mapfile)
	{
		if(isBoundaryMaps == false) {
			cout << "UMesh2d: writeBoundaryMapsToFile(): ! Boundary maps not available!" << endl;
			return;
		}
		ofstream ofile(mapfile);
		ofile << nbface << '\n'<< "bifmap\n";
		for(int i = 0; i < nbface; i++)
			ofile << bifmap.get(i) << ' ';
		ofile << '\n';
		ofile << "ifbmap\n";
		for(int i = 0; i < nbface; i++)
			ofile << ifbmap.get(i) << ' ';
		ofile << '\n';
		ofile.close();
	}

	void readBoundaryMapsFromFile(string mapfile)
	{
		ifstream ofile(mapfile);
		string dum; int sz;
		ofile >> sz >>  dum;
		cout << "UMesh2d: readBoundaryMapsFromFile(): Number of boundary faces in file = " << sz << endl;
		bifmap.setup(sz,1);
		ifbmap.setup(sz,1);

		for(int i = 0; i < nbface; i++)
			ofile >> bifmap(i);

		ofile >> dum;
		for(int i = 0; i < nbface; i++)
			ofile >> ifbmap(i);

		ofile.close();
		isBoundaryMaps = true;
	}

	void compute_intfacbtags()
	{
		/// Populate intfacbtags with boundary markers of corresponding bfaces

		intfacbtags.setup(nbpoin,nbtag);

		if(isBoundaryMaps == false)
		{
			cout << "UMesh2d: compute_intfacbtags(): ! Boundary maps are not available!" << endl;
			return;
		}

		for(int ibface = 0; ibface < nface; ibface++)
		{
			for(int j = 0; j < nbtag; j++)
				intfacbtags(ifbmap(ibface),j) = bface(ibface,nnofa+j);
		}
	}

	UMesh2d convertLinearToQuadratic()
	{
		cout << "UMesh2d: convertLinearToQuadratic(): Producing quadratic mesh from linear mesh" << endl;
		UMesh2d q;
		if(nnofa != 2) { cout << "! UMesh2d: convertLinearToQuadratic(): Mesh is not linear!!" << endl; return q;}

		if(nnode == 3)			// for simplicial mesh
		{
			cout << "UMesh2d: convertLinearToQuadratic(): Simplicial mesh." << endl;
			int parm = 1;		// 1 extra node per face
			q.ndim = ndim;
			q.npoin = npoin + naface;
			q.nelem = nelem;
			q.nface = nface;
			q.nbface = nbface;
			q.naface = naface;
			q.nnofa = nnofa+parm;
			q.nnode = nnode + nfael*parm;
			q.nfael = nfael;
			q.nbtag = nbtag;
			q.ndtag = ndtag;

			q.coords.setup(q.npoin, q.ndim);
			q.inpoel.setup(q.nelem, q.nnode);
			q.bface.setup(q.nface, q.nnofa+q.nbtag);

			for(int i = 0; i < npoin; i++)
				for(int j = 0; j < ndim; j++)
					q.coords(i,j) = coords(i,j);

			for(int i = 0; i < nelem; i++)
				for(int j = 0; j < nnode; j++)
					q.inpoel(i,j) = inpoel(i,j);

			for(int i = 0; i < nface; i++)
			{
				for(int j = 0; j < nnofa; j++)
					q.bface(i,j) = bface(i,j);
				for(int j = nnofa; j < nnofa+nbtag; j++)
					q.bface(i,j+parm) = bface(i,j);
			}

			q.vol_regions = vol_regions;

			int ied, p1, p2, ielem, jelem, idim, inode, lp1, lp2, ifa;

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over boundary faces..." << endl;
			// iterate over boundary faces
			for(ied = 0; ied < nbface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+ied;

				// find the bface that this face corresponds to
				for(ifa = 0; ifa < nface; ifa++)
				{
					if((p1 == bface(ifa,0) && p2 == bface(ifa,1)) || (p1 == bface(ifa,1) && p2 == bface(ifa,0)))	// face found
					{
						q.bface(ifa,nnofa) = npoin+ied;
					}
				}
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over internal faces..." << endl;
			// iterate over internal faces
			for(ied = nbface; ied < naface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+ied;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(jelem,inode)) lp1 = inode;
					if(p2 == inpoel(jelem,inode)) lp2 = inode;
				}

				// in the right element, the new point is in face ip2
				q.inpoel(jelem, nnode+lp2) = npoin+ied;
			}
			cout << "UMesh2d: convertLinearToQuadratic(): Done." << endl;
			return q;
		}
		else 					// for non-simplicial mesh, add extra points at cell-centres as well
		{
			cout << "UMesh2d: convertLinearToQuadratic(): Non-simplicial mesh." << endl;
			int parm = 1;		// 1 extra node per face
			q.ndim = ndim;
			q.npoin = npoin + naface + nelem;
			q.nelem = nelem;
			q.nface = nface;
			q.nbface = nbface;
			q.naface = naface;
			q.nnofa = nnofa+parm;
			q.nnode = nnode + nfael*parm + 1;
			q.nfael = nfael;
			q.nbtag = nbtag;
			q.ndtag = ndtag;

			q.coords.setup(q.npoin, q.ndim);
			q.inpoel.setup(q.nelem, q.nnode);
			q.bface.setup(q.nface, q.nnofa+q.nbtag);

			for(int i = 0; i < npoin; i++)
				for(int j = 0; j < ndim; j++)
					q.coords(i,j) = coords(i,j);

			for(int i = 0; i < nelem; i++)
				for(int j = 0; j < nnode; j++)
					q.inpoel(i,j) = inpoel(i,j);

			for(int i = 0; i < nface; i++)
			{
				for(int j = 0; j < nnofa; j++)
					q.bface(i,j) = bface(i,j);
				for(int j = nnofa; j < nnofa+nbtag; j++)
					q.bface(i,j+parm) = bface(i,j);
			}

			q.vol_regions = vol_regions;

			int ied, p1, p2, ielem, jelem, idim, inode, lp1, lp2, ifa;

			// get cell centres
			for(int iel = 0; iel < nelem; iel++)
			{
				double c_x = 0, c_y = 0;

				for(int inode = 0; inode < nnode; inode++)
				{
					c_x += coords(inpoel(iel,inode),0);
					c_y += coords(inpoel(iel,inode),1);
				}
				c_x /= nnode;
				c_y /= nnode;
				q.coords(npoin+iel,0) = c_x;
				q.coords(npoin+iel,1) = c_y;
				q.inpoel(iel,q.nnode-1) = npoin+iel;
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over boundary faces..." << endl;
			// iterate over boundary faces
			for(ied = 0; ied < nbface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+nelem+ied;

				// find the bface that this face corresponds to
				for(ifa = 0; ifa < nface; ifa++)
				{
					if((p1 == bface(ifa,0) && p2 == bface(ifa,1)) || (p1 == bface(ifa,1) && p2 == bface(ifa,0)))	// face found
					{
						q.bface(ifa,nnofa) = npoin+nelem+ied;
					}
				}
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over internal faces..." << endl;
			// iterate over internal faces
			for(ied = nbface; ied < naface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+nelem+ied;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(jelem,inode)) lp1 = inode;
					if(p2 == inpoel(jelem,inode)) lp2 = inode;
				}

				// in the right element, the new point is in face ip2
				q.inpoel(jelem, nnode+lp2) = npoin+nelem+ied;
			}
			cout << "UMesh2d: convertLinearToQuadratic(): Done." << endl;
			return q;
		}
	}
};

} // end namespace acfd