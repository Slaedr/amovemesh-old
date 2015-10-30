/** Class to convert linear mesh into quadratic mesh using P2 linear elasticity.
	Aditya Kashi
	October 29, 2015
*/

#ifndef __AGEOMETRY_H
	#include <ageometry.hpp>
#endif

#ifndef __ALINALG_H
	#include <alinalg.hpp>
#endif

#ifndef __ALINELAST_P2_SPARSE_H
	#include <alinelast_p2_sparse.hpp>
#endif

#define __ALINELAST_CURVEDMESHGEN2D_H 1

using namespace std;
using namespace amat;
using namespace acfd;


/** Class to generate curved mesh from a linear mesh using cubic spline reconstruction and one of the mesh movement techniques. */

class Curvedmeshgen2d
{
	UMesh2d* m;						///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh2d* mq;					///< Data of the corresponding (straight-faced) quadratic mesh
	LinElastP2* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction2d br;	///< Object to reconstruct the boundary using cubic splines.
	
	double tol;						///< Tolerance for linear solvers.
	double maxiter;					///< Maximum number of iterations for linear solvers.
	double young;					///< Young's modulus
	double nu;						///< Poisson's ratio
	double lambda;
	double mu;

	int nbounpoin;					///< Number if boundary points.
	int ninpoin;					///< Number of interior points.
	Matrix<double> disps;			///< Displacement of midpoint of each face
	Matrix<double> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	Matrix<double> bounpoints;
	Matrix<double> inpoints;
	Matrix<int> bflagg;				///< This flag is true if the corresponding mesh node lies on a boundary.
	Matrix<int> toRec;				///< This flag is true if a boundary face is to be reconstructed.

public:
	void setup(UMesh2d* mesh, UMesh2d* meshq, LinElastP2* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double toler, double maxitera, double youngsmodulus, double poissonsratio); 

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void Curvedmeshgen2d::setup(UMesh2d* mesh, UMesh2d* meshq, LinElastP2* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double toler, double maxitera, double youngsmodulus, double poissonsratio)
{
	m = mesh;
	mq = meshq;
	mmv = mmove;
	br.setup(m, num_parts, boundarymarkers, angle_threshold);
	tol = toler;
	maxiter = maxitera;
	young = youngsmodulus;
	nu = poissonsratio;
	lambda = nu*young/((1+nu)*(1-2*nu));
	mu = young/(2*(1+nu));

	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);
	
	// demarcate which faces are to be reconstructed
	toRec.setup(m->gnface(),1);
	toRec.zeros();
	for(int iface = 0; iface < m->gnface(); iface++)
		for(int i = 0; i < boundarymarkers.size(); i++)
			for(int j = 0; j < boundarymarkers[i].size(); j++)
				if(m->gbface(iface,m->gnnofa()) == boundarymarkers[i][j])
					toRec(iface) = 1;
}

/** Computes displacement of midpoint of each face. */
void Curvedmeshgen2d::compute_boundary_displacements()
{
	br.preprocess();
	br.detect_corners();
	br.split_parts();
	br.compute_splines(tol,maxiter);

	// get coords of midpoints of each face
	Matrix<double> facemidpoints(m->gnface(),m->gndim());
	for(int iface = 0; iface < m->gnface(); iface++)
		for(int idim = 0; idim < m->gndim(); idim++)
		{
			double sum = 0;
			for(int inode = 0; inode < m->gnnofa(); inode++)
				sum += m->gcoords(m->gbface(iface,inode),idim);
			facemidpoints(iface,idim) = sum/m->gnnofa();
		}
	
	double uh = 0.5;
	for(int iface = 0; iface < m->gnface(); iface++)
	{
		// first check if iface was reconstructed!
		if(toRec(iface))
				for(int idim = 0; idim < m->gndim(); idim++)
					disps(iface,idim) = br.getcoords(iface,idim,uh) - facemidpoints.get(iface,idim);
	}

	/// We do not need the linear mesh once we have the displacements of the faces' midpoints.
}

/** Uses the previously computed displacements of the face midpoints to curve the mesh.
*/
void Curvedmeshgen2d::generate_curved_mesh()
{
	/** 
	Note that this function works with the straight quadratic mesh.
	We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
	*/

	mmv->setup(mq, mu, lambda);
	cout << "Cuvedmeshgen2d: generate_curved_mesh(): Assembling stiffness matrix and load vector." << endl;
	mmv->assembleStiffnessMatrix();
	mmv->assembleLoadVector();
	
	cout << "Cuvedmeshgen2d: generate_curved_mesh(): Applying Dirichlet BCs." << endl;
	mmv->dirichletBC_onAllBface(disps);

	Matrix<double> alldisps(mq->gnpoin()*2,1);

	SpMatrix A = mmv->stiffnessMatrix();
	Matrix<double> b = mmv->loadVector();

	Matrix<double> xold(mq->gnpoin(),1);
	xold.zeros();

	alldisps = sparseCG_d(&A, b, xold, tol, maxiter);

	Matrix<double> fpos(mq->gnpoin(),mq->gndim());
	for(int i = 0; i < mq->gnpoin(); i++)
		for(int j = 0; j < mq->gndim(); j++)
			fpos(i,j) = mq->gcoords(i,j);

	// update poistions
	for(int j = 0; j < mq->gndim(); j++)
	{
		for(int i = 0; i < mq->gnpoin(); i++)
			fpos(i,j) += alldisps(j*mq->gnpoin() + i);
	}
	//update mesh
	mq->setcoords(&fpos);
	
}

// ------------ end --------------------
