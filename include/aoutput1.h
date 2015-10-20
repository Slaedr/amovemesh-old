#ifndef __AMATRIX2_H
#include "amatrix2.h"
#endif

#ifndef __AMESH_H
#include "amesh2.h"
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
		out << "\t\t\t\t" << m.ginpoel(i,0)-1 << " " << m.ginpoel(i,1)-1 << " " << m.ginpoel(i,2)-1 << '\n';
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
		out << "\t\t\t\t" << m.ginpoel(i,0)-1 << " " << m.ginpoel(i,1)-1 << " " << m.ginpoel(i,2)-1 << '\n';
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
		out << "\t\t\t\t" << m.ginpoel(i,0)-1 << " " << m.ginpoel(i,1)-1 << " " << m.ginpoel(i,2)-1 << '\n';
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

// currently works for 2-d meshes only
void writeMeshGmsh2(UTriMesh m, string fname)
{
	ofstream out(fname.c_str());
	out << "$MeshFormat \n2.2 0 8 \n$EndMeshFormat\n";
	out << "$Nodes\n"; 
	out << m.gnnode() << '\n';
	for(int i = 0; i < m.gnnode(); i++)
	{
		out << i+1;
		for(int j = 0; j < m.gndim(); j++)
			out << " " << m.gcoords(i,j);
		out << " " << 0 << '\n';
	}
	out << "$EndNodes\n";
	
	out << "$Elements\n";
	out << m.gnelem()+m.gnface() << '\n';
	int iel = 0;
	
	out.close();
}
