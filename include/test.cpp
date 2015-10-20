#include "asparsematrix.hpp"
#include "alinalg.hpp"

using namespace std;
using namespace amat;

int main()
{
	MatrixCOO A(6,6);
	Matrix<double> prod(6,1);
	//cout << A.get(1,1) << endl;
	//A.set(1,1,10.5);
	A.set(0,0, 5.1);
	A.set(0,1, 2);
	A.set(1,0, 2.05);
	A.set(1,1, 4);
	A.set(2,2, 9);
	A.set(2,3, -1.5);
	A.set(3,2, -1.52);
	A.set(3,5, 2);
	A.set(5,3, 2.1);
	A.set(3,3, 6);
	A.set(4,4, 2.5);
	A.set(4,1, -0.5);
	A.set(1,4, -0.48);
	A.set(5,5, 8);
	A.set(5,2, -2);
	A.set(2,5, -2.005);
	A.set(5,4, 1.5);
	A.set(4,5, 1.55);

	A.mprint();

	Matrix<double> x(6,1);
	x(0) = 9; x(1) = 8; x(2) = 7; x(3) = 6; x(4) = 5; x(5) = 4;
	x.mprint();
	Matrix<double> y(6,1);
	y.ones();
	cout << "Data printed." << endl;

	Matrix<double> xold(6,1);
	xold.zeros();

	//Matrix<double> ans = sparsegaussseidel(&A, x, xold, 1e-6, 600);
	//Matrix<double> ans = sparseSOR(&A, x, xold, 1e-6, 600, 1.25);
	//Matrix<double> ans = sparseCG_d(&A, x, xold, 1e-6, 600);
	Matrix<double> ans = sparse_bicgstab(&A, x, xold, 1e-6, 600);
	ans.mprint();

	A.multiply(ans, &prod);

	cout << "error is " << (prod - x).l2norm() << endl;

	prod.mprint();

	/*cout << "Matrix multiply\n";
	A.multiply(x, &prod);
	prod.mprint();*/

	cout << endl;
	return 0;
}
