
#include "../include/UsefulFunctions.h"

/**************************************************************************/
/*  UsefulFunctions.cc                                           		  */
/*  Libreria in cui definisco funzioni ausiliarie						  */
/**************************************************************************/

void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_point_data(meshFEM, solution, solutionName.data());

    return;

}// exportSolution


void exportMesh ( const std::string& fileName, const getfem::mesh& mesh )
{
    getfem::vtk_export exporter(fileName.data(), true);

    // con true mi permette di esportare la mesh in formato txt

    exporter.exporting(mesh);

    exporter.write_mesh();

    return;
}// exportMesh


scalar_type min ( const std::string& flusso, const scalar_type ul, const scalar_type ur )
{
	LifeV::Parser M_parser;
	scalar_type tmp;
	scalarVector_Type u( 500 );
	scalar_type h = ( ur-ul )/500;

	for ( size_type i = 0; i < u.size(); i++ )
	{
		u[ i ] = ul + h*i;
	}

	M_parser.setString ( flusso );
	M_parser.setVariable ( "x", u [ 0 ] );

	tmp=M_parser.evaluate();

	for( size_type i = 0; i < u.size(); i++ )
	{
		M_parser.setVariable ( "x", u [ i ] );

		if ( M_parser.evaluate() < tmp )
		{
			tmp = M_parser.evaluate();
		}
	}

	return tmp;
}// min



scalar_type max ( const std::string& flusso, const scalar_type ul, const scalar_type ur )
{
	LifeV::Parser M_parser;
	scalar_type tmp;
	scalarVector_Type u( 500 );
	scalar_type h = ( ul-ur )/500;

	for ( size_type i = 0; i < u.size(); i++ )
	{
		u[ i ] = ur + h*i;
	}

	M_parser.setString ( flusso );
	M_parser.setVariable ( "x", u [ 0 ] );

	tmp=M_parser.evaluate();

	for( size_type i = 0; i < u.size(); i++ )
	{
		M_parser.setVariable ( "x", u [ i ] );

		if ( M_parser.evaluate() > tmp )
		{
			tmp = M_parser.evaluate();
		}
	}

	return tmp;

}// max


VectorIntersection::VectorIntersection ()
{}


