#include <amesh2d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string dum, mfile, omfile;
	
	mfile = "feflo.domn.cylinder.coarse";
	omfile = "cylinder-coarse.msh";

	UMesh2d m;
	m.readDomn(mfile);
	m.writeGmsh2(omfile);

	cout << endl;
	return 0;
}
