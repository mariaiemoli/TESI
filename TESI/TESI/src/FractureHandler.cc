
#include "../include/FractureHandler.h"

/**************************************************************************/
/*  FractureHandler.cc                                                	  */
/*  Classe che inizializza e gestisce una frattura				 		  */
/**************************************************************************/

FractureHandler::FractureHandler ( const GetPot& dataFile,
								   const size_type ID,
								   const ExporterPtr_Type& exporter,
								   const std::string& section,
								   const std::string& sectionSaturation):
								   M_ID( ID ),
								   M_exporter ( exporter ),
								   M_data( dataFile, section )
{
	M_levelSet.reset( new LevelSetHandler_Type ( dataFile, section ) );

}// costruttore



void FractureHandler::init ()
{
    // Geometric transformation usign primal finite elements type in the fracture
    M_geometricTransformation = bgeot::geometric_trans_descriptor( M_data.getMeshType() );


    //-------------------- M_meshFlat, mesh sul piano orizzontale --------------------//

    // costruisco una mesh standard, l'idea di mantenere una mesh 1d sul piano orizzontale è solo per comodità

    sizeVector_Type fractureNumberSubdivision( M_data.getSpaceDimension() );	// numero suddivisioni

 //   std::cout << "  M_data.getSpatialDiscretization():  " <<  M_data.getSpatialDiscretization() << std::endl;

    std::fill( fractureNumberSubdivision.begin(), fractureNumberSubdivision.end(), M_data.getSpatialDiscretization() );

    getfem::regular_unit_mesh( M_meshFlat, fractureNumberSubdivision, M_geometricTransformation );

    bgeot::base_matrix fractureTransformMatrix( M_data.getSpaceDimension(), M_data.getSpaceDimension() );

    scalarVector_Type fractureLength(3);
    fractureLength [ 0 ] = M_data.getLengthAbscissa();
    fractureLength [ 1 ] = M_data.getLengthOrdinate();
    fractureLength [ 2 ] = M_data.getLengthQuota();

    for ( size_type i = 0; i < M_data.getSpaceDimension(); ++i )
    {
        fractureTransformMatrix(i, i) = (i < M_data.getSpaceDimension()) ? fractureLength [ i ] : 1.0;
    }

    // riscalo la mesh unitaria M_mediumMesh to [M_lengthAbscissa,M_lengthOrdinate]
    M_meshFlat.transformation( fractureTransformMatrix );
    for ( size_type i = 0; i < M_meshFlat.points().size(); ++i )
    {
        M_meshFlat.points() [ i ] [ 0 ] = M_data.meshSpacing( M_meshFlat.points() [ i ] [ 0 ] );
    }

    /*
    std::ostringstream osFileName;

    osFileName << "cmeshFlat-" << ".vtk";
    exportMesh(M_exporter->getFolder() + osFileName.str().c_str(), M_meshFlat );
	*/


    //-------------------- M_mesh, mesh reale nel piano 2d --------------------//
    //costruisco la M_mediumMesh n.2 quella di servizio, quella della frattura mappata



    sizeVector_Type ind( M_meshFlat.points().size() );
    for ( size_type i = 0; i < M_meshFlat.points().size(); ++i )
    {
        bgeot::base_node P(M_data.getSpaceDimension() + 2);
        bgeot::base_node P1(M_data.getSpaceDimension() + 2);

        scalar_type t = i*1./(M_data.getSpatialDiscretization () );

        P1 [ 0 ] = t;

        P [ 0 ] = P1 [ 0 ];
        P [ 1 ] =  M_levelSet->getData()->y_map( P1 );
        ind [ i ] = M_mesh.add_point( P );
    }



    for ( size_type i = 0; i < M_meshFlat.convex_index().size(); ++i )
    {
        std::vector<bgeot::size_type> point(3);
        point [ 0 ] = M_meshFlat.ind_points_of_convex(i) [ 0 ];
        point [ 1 ] = M_meshFlat.ind_points_of_convex(i) [ 1 ];
        point [ 2 ] = M_meshFlat.ind_points_of_convex(i) [ 2 ];

        M_mesh.add_convex( M_geometricTransformation, point.begin() );
    }

    /*
    std::cout << " esporto la seconda mesh " << std::endl;
    osFileName << "cmesh-" << ".vtk";
	exportMesh(M_exporter->getFolder() + osFileName.str().c_str(), M_mesh );

	std::cout << " fatto " << std::endl;
 	 */


    // inizializzo M_ci, condizione iniziale
	M_ci.resize( M_data.getSpatialDiscretization() );

//	std::cout << " M_data.getX0():  " << M_data.getX0() << std::endl;

	for ( size_type i = 0; i < M_data.getSpatialDiscretization(); i++ )
	{
		bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_meshFlat.points_of_convex ( i );

		scalar_type x = nodes [ 0 ] [ 0 ];

		if ( x <= M_data.getX0() )
		{
			M_ci [ i ] = M_data.getUl();
		}
		else
		{
			M_ci [ i ] = M_data.getUr();
		}
	}


	std::ofstream exp("u0.txt");

			exp << M_ci;


	return;
}// init









