
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
								   M_data( dataFile, section ),
								   M_exporter ( exporter )
{
	M_levelSet.reset( new LevelSetHandler_Type ( dataFile, section ) );

}// costruttore



void FractureHandler::init ()
{
    // Geometric transformation usign primal finite elements type in the fracture
    M_geometricTransformation = bgeot::geometric_trans_descriptor( M_data.getMeshType() );


    //-------------------- M_meshFlat, mesh sul piano orizzontale --------------------//

    // costruisco una mesh standard, l'idea di mantenere una mesh 1d sul piano orizzontale è solo per comodità

    sizeVector_Type fractureNumberSubdivision( M_data.getSpaceDimension() );

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

    // inizializzo M_ci, condizione iniziale
	M_ci.resize( M_data.getSpatialDiscretization() - 1 );



	for ( size_type i = 0; i < M_data.getSpatialDiscretization() - 1; i++ )
	{
		bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_meshFlat.points_of_convex ( i );

		scalar_type x = nodes [ 0 ] [ 0 ];

		M_ci [ i ] = M_data.getCi( x );
	}

	std::ostringstream ss;
	ss << "./matlab/risultati/inizialSaturation_" << M_ID << ".txt";

	std::string name = ss.str();

	std::ofstream exp( name );

	exp << M_ci;

	return;
}// init


void FractureHandler::setFractureIntersection ( const sizeVector_Type& nodes, const FractureHandler& fractureInvolved )
{
	size_type otherFractureID = fractureInvolved.getID();

	M_fractureIntersectElements[ otherFractureID ].push_back( nodes[ 0 ] );

	return;
}// setFractureIntersection

