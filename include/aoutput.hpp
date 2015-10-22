/** A collection of subroutines to write mesh data to various kinds of output formats.*/

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#ifndef __AMESH2D_H
#include "amesh2.hpp"
#endif

#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
#endif

#ifndef __AMESH2DCURVED_H
#include "amesh_curved.hpp"
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_STRING
#include <string>
#endif

#define __AOUTPUT_H 1

using namespace std;
using namespace amat;
using namespace acfd;

void writeScalarsVectorToVtu_CellData(string fname, const UMesh2d& m, const Matrix<double>& x, string scaname[], const Matrix<double>& y, string vecname);

// -------------- Implemetations -----------------------

void prepareVtk (ofstream& out, UTriMesh m, string title)
{
	string start = "# vtk DataFile Version 2.0";
	out << start << '\n';
	out << title << '\n';
	out << "ASCII" << '\n';
	out << '\n';

	start = "DATASET UNSTRUCTURED_GRID";
	out << start << '\n';

	out << "POINTS " << m.gnpoin() << " double" << '\n';
	for(int i = 0; i < m.gnpoin(); i++)
		out << "  " << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << '\n';

	int cellrowsize = m.gnnode()*m.gnelem();

	out << "CELLS " << m.gnelem() << " " << cellrowsize << '\n';
	for(int i = 0; i < m.gnelem(); i++)
	{
		out << "  " << m.gnnode();
		for(int j = 0; j < m.gnnode(); j++)
			out << ' ' << m.ginpoel(i,j);
		out << '\n';
	}

	// In vtk, linear interpolation cell is type 5. Quadratic interpolation cell is type 22
	int celltype = 5;
	out << '\n';
	out << "CELL_TYPES " << m.gnelem() << '\n';
	for(int i = 0; i < m.gnelem(); i++)
		out << celltype << '\n';
	out << '\n';


}

void writeScalarToVtk(ofstream& out, UTriMesh m, Matrix<double> x, string scaname)
{
	out << "POINT_DATA " << m.gnpoin() << '\n';
	out << "SCALARS " << scaname << " double" << '\n';
	out << "LOOKUP_TABLE default\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "  " << x(i,0) << '\n';
		//out << "  " << i << '\n';
	out << "\n";
}

void writeScalarToVtu(string fname, UTriMesh m, Matrix<double> x, string scaname)
{
	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter point scalar data
	out << "\t\t<PointData Scalars=\"scalars\">\n";
	out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname << "\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t\t" << x(i,0) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</PointData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.ginpoel(i,0) << " " << m.ginpoel(i,1) << " " << m.ginpoel(i,2) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 5 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

// writes a single scalar and vector
void writeScalarVectorToVtu(string fname, UTriMesh m, Matrix<double> x, string scaname, Matrix<double> y, string vecname)
{
	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter point scalar data
	out << "\t\t<PointData Scalars=\""<<scaname<< "\" Vectors=\"" << vecname << "\">\n";
	out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname << "\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t\t" << x(i,0) << '\n';
	out << "\t\t\t</DataArray>\n";
	//out << "\t\t</PointData>\n";

	//enter vector point data
	//out << "\t\t<PointData Vectors=\"vectors\">\n";
	out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << vecname << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t\t" << y(i,0) << " " << y(i,1) << " " << 0.0 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</PointData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.ginpoel(i,0) << " " << m.ginpoel(i,1) << " " << m.ginpoel(i,2) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 5 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

void writeScalarsVectorToVtu_CellData(string fname, UTriMesh& m, Matrix<double>& x, string scaname[], Matrix<double>& y, string vecname)
{
	cout << "aoutput: Writing vtu output...\n";
	ofstream out(fname);

	const int celltype = 5;

	int nscalars = x.cols();

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter point scalar data
	out << "\t\t<CellData Scalars=\""<<scaname[0]<< "\" Vectors=\"" << vecname << "\">\n";

	//cout << "aoutput: Writing scalars..\n";
	for(int in = 0; in < nscalars; in++)
	{
		out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname[in] << "\" Format=\"ascii\">\n";
		for(int i = 0; i < m.gnelem(); i++)
			out << "\t\t\t\t" << x(i,in) << '\n';
		out << "\t\t\t</DataArray>\n";
	}
	//cout << "aoutput: Scalars written.\n";

	//enter vector point data
	out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << vecname << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << y(i,0) << " " << y(i,1) << " " << 0.0 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</CellData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.ginpoel(i,0) << " " << m.ginpoel(i,1) << " " << m.ginpoel(i,2) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << celltype << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

/** Writes multiple scalar data sets and one vector data set, all cell-centered data, to a file in VTU format.
	If either x or y is a 0x0 matrix, it is ignored.
*/
void writeScalarsVectorToVtu_CellData(string fname, const UMesh2d& m, const Matrix<double>& x, string scaname[], const Matrix<double>& y, string vecname)
{
	int elemcode = 5;
	if(m.gnnode() == 4)
		elemcode = 9;
	else if(m.gnnode() == 6)
		elemcode = 22;
	else if(m.gnnode() == 8)
		elemcode = 23;
	else if(m.gnnode() == 9)
		elemcode = 28;
	
	cout << "aoutput: Writing vtu output...\n";
	ofstream out(fname);

	int nscalars = x.cols();

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//out << "\t\t<CellData Scalars=\""<<scaname[0]<< "\" Vectors=\"" << vecname << "\">\n";
	
	if(x.msize()>0 || y.msize()>0) {
		out << "\t\t<CellData ";
		if(x.msize() > 0)
			out << "Scalars=\"" << scaname[0] << "\" ";
		if(y.msize() > 0)
			out << "Vectors=\"" << vecname << "\"";
		out << ">\n";
	}
	
	//enter cell scalar data if available
	if(x.msize() > 0) {
		//cout << "aoutput: Writing scalars..\n";
		for(int in = 0; in < nscalars; in++)
		{
			out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname[in] << "\" Format=\"ascii\">\n";
			for(int i = 0; i < m.gnelem(); i++)
				out << "\t\t\t\t" << x.get(i,in) << '\n';
			out << "\t\t\t</DataArray>\n";
		}
		//cout << "aoutput: Scalars written.\n";
	}

	//enter vector cell data if available
	if(y.msize() > 0) {
		out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << vecname << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		for(int i = 0; i < m.gnelem(); i++)
		{
			out << "\t\t\t\t";
			for(int idim = 0; idim < y.cols(); idim++)
				out << y.get(i,idim) << " ";
			if(y.cols() == 2)
				out << "0.0 ";
			cout << '\n';
		}
		out << "\t\t\t</DataArray>\n";
	}
	if(x.msize() > 0 || y.msize() > 0)
		out << "\t\t</CellData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
	{
		out << "\t\t\t";
		for(int idim = 0; idim < m.gndim(); idim++)
			out << m.gcoords(i,idim) << " ";
		if(m.gndim() == 2)
			out << "0.0 ";
		out << '\n';
	}
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++) {
		out << "\t\t\t\t"; 
		for(int inode = 0; inode < m.gnnode(); inode++)	
			out << m.ginpoel(i,inode) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.gnnode()*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << elemcode << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

void writeScalarsVectorToVtu(string fname, UTriMesh& m, Matrix<double> x, string scaname[], Matrix<double> y, string vecname)
{
	cout << "aoutput: Writing vtu output...\n";
	ofstream out(fname);

	const int celltype = 5;

	int nscalars = x.cols();

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter point scalar data
	out << "\t\t<PointData Scalars=\""<<scaname[0]<< "\" Vectors=\"" << vecname << "\">\n";

	//cout << "aoutput: Writing scalars..\n";
	for(int in = 0; in < nscalars; in++)
	{
		out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname[in] << "\" Format=\"ascii\">\n";
		for(int i = 0; i < m.gnpoin(); i++)
			out << "\t\t\t\t" << x(i,in) << '\n';
		out << "\t\t\t</DataArray>\n";
	}
	//cout << "aoutput: Scalars written.\n";

	//enter vector point data
	out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << vecname << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t\t" << y(i,0) << " " << y(i,1) << " " << 0.0 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</PointData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.ginpoel(i,0) << " " << m.ginpoel(i,1) << " " << m.ginpoel(i,2) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << celltype << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

void writecsv (string fname, UTriMesh m, Matrix<double> x, string scaname)
{
	ofstream out(fname.c_str());

	out << "x-coords, y-coords, " << scaname << '\n';

	for(int i = 0; i < m.gnpoin(); i++)
		out << "  " << m.gcoords(i,0) << ", " << m.gcoords(i,1) << ", " << x(i,0) << '\n';

	out.close();
}

void writeMeshToVtu(string fname, UTriMesh m)
{
	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.ginpoel(i,0) << " " << m.ginpoel(i,1) << " " << m.ginpoel(i,2) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 5 << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

void writeQuadraticMeshToVtu(string fname, UTriMeshCurved& m)
{
	int elemcode = 22;
	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
	{
		out << "\t\t\t\t";
		for(int j = 0; j < m.gnnode(); j++)
			out << m.ginpoel(i,j) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << 3*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << elemcode << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

void writeQuadraticMeshToVtu(string fname, UMesh2d& m)
{
	int elemcode = 5;
	if(m.gnnode() == 4)
		elemcode = 9;
	else if(m.gnnode() == 6)
		elemcode = 22;
	else if(m.gnnode() == 8)
		elemcode = 23;
	else if(m.gnnode() == 9)
		elemcode = 23;

	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
	{
		out << "\t\t\t\t";
		for(int j = 0; j < m.gnnode(); j++)
			out << m.ginpoel(i,j) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.gnnode()*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << elemcode << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}
