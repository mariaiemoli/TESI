
#include "../include/FractureHandler.h"

/**************************************************************************/
/*  FractureHandler.cc                                                	  */
/*  Classe che inizializza e gestisce una frattura				 		  */
/**************************************************************************/

FractureHandler::FractureHandler ( const GetPot& dataFile,
								   const size_type ID,
								   const ExporterPtr_Type& exporter,
								   const std::string& section ):
								   M_ID( ID ),
								   M_data( dataFile, section ),
								   M_exporter ( exporter ),
                                   M_meshFEM( M_meshFlat ),
                                   M_meshFEM2( M_meshFlat ),
                                   M_meshFEMVisualization( M_mesh ),
                                   M_integrationMethod( M_meshFlat ),
                                   M_integrationMethodVisualization( M_mesh ),
                                   M_integrationMethod2( M_meshFlat ),
                                   M_integrationMethodLinear( M_mesh ),
                                   M_meshFEMLinear( M_mesh )

{
	M_levelSet.reset( new LevelSetHandler_Type ( dataFile, section ) );

	base_node t0(1);
	base_node t1(1);
	t0[ 0 ] = 0;
	t1[ 0 ] = 1;

	scalar_type ax = M_data.getA();
	scalar_type ay = M_levelSet->getData ( )->y_map( t0 );
	scalar_type bx = M_data.getB();
	scalar_type by = M_levelSet->getData ( )->y_map( t1 );

	scalar_type spatialDiscretization = M_data.getSpatialDiscretization ();

	M_h = sqrt( (ax-bx)*(ax-bx) + (ay-by)*(ay-by) )/( spatialDiscretization-1 );

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

        P [ 0 ] = M_meshFlat.points() [ i ] [ 0 ];

        scalar_type t = i*1./(M_data.getSpatialDiscretization () );

        P1 [ 0 ] = t;
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

    // Definisco gli elementi finiti

    // Definisco gli elementi finiti per la velocità nella frattura
    getfem::pfem fractureFEType2 = getfem::fem_descriptor( M_data.getFEMType2() );

    // Definisco il tipo di integrazione per la velocità nella frattura
    getfem::pintegration_method fractureIntegrationType2 = getfem::int_method_descriptor( M_data.getIntegrationType2() );

    // Definisco il metodo di integrazione per la velocità nella frattura
    M_integrationMethod2.set_integration_method( M_meshFlat.convex_index(), fractureIntegrationType2 );

    // Definisco lo spazio degli elementi finiti per la velocità nella frattura
    M_meshFEM2.set_qdim( M_data.getSpaceDimension() );
    M_meshFEM2.set_finite_element( M_meshFlat.convex_index(), fractureFEType2 );

    // Definisco gli elementi finiti per i coefficienti nella frattura
    getfem::pfem fractureFEMTypeLinear = getfem::fem_descriptor( M_data.getFEMTypeLinear() );

    // Definisco il tipo di integrazione per i coefficienti nella frattura
    getfem::pintegration_method fractureIntegrationTypeLinea = getfem::int_method_descriptor( M_data.getIntegrationType2() );

    // Definisco il metodo di integrazione per i coefficienti nella frattura
    M_integrationMethodLinear.set_integration_method( M_mesh.convex_index(), fractureIntegrationTypeLinea );

    // Definisco lo spazio degli elementi finiti per i coefficienti nella frattura
    M_meshFEMLinear.set_qdim( M_data.getSpaceDimension() );
    M_meshFEMLinear.set_finite_element( M_mesh.convex_index(), fractureFEMTypeLinear );

    // P(K-1) discontinuous FE
    // Definisco gli elementi finiti nella frattura
    getfem::pfem fractureFEType = getfem::fem_descriptor( M_data.getFEMType () );

    // Definisco il tipo di integrazione nella frattura
    getfem::pintegration_method fractureIntegrationType = getfem::int_method_descriptor( M_data.getIntegrationType() );

    // Definisco il metodo di integrazione nella frattura
    M_integrationMethod.set_integration_method( M_meshFlat.convex_index(), fractureIntegrationType );

    // Definisco lo spazio degli elementi finiti per i coefficienti nella frattura
    M_meshFEM.set_finite_element(M_meshFlat.convex_index(), fractureFEType );

    // Metodo di integrazione per la pressione nella frattura per la visualizzazione
    M_integrationMethodVisualization.set_integration_method( M_mesh.convex_index(), fractureIntegrationType );

    // Elementi finiti per la pressione nella frattura per la visualizzazione
    M_meshFEMVisualization.set_finite_element( M_mesh.convex_index(), fractureFEType );

    // Allocate data

    // Alloco il vettore per M_etaNormalInterpolated
    gmm::resize ( M_etaNormalInterpolated, M_meshFEM.nb_dof() );
    gmm::clear ( M_etaNormalInterpolated );

    // Alloco il vettore per M_etaTangentialInterpolated
    gmm::resize ( M_etaTangentialInterpolated, M_meshFEM.nb_dof() );
    gmm::clear ( M_etaTangentialInterpolated );

    // Riempio i vettori M_etaNormalInterpolated, M_etaTangentialInterpolated, M_muNormalInterpolated e M_muTangentialInterpolated della frattura
    for ( size_type i = 0; i < M_meshFEM.nb_dof(); ++i )
    {
        M_etaNormalInterpolated [ i ] = M_data.etaNormalDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaNormal();

        M_etaTangentialInterpolated [ i ] = M_data.etaTangentialDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaTangential();
    }

    gmm::resize ( M_inverseMeshSize, M_meshFEM.nb_dof() );
    gmm::clear ( M_inverseMeshSize );

    M_meshFlat.region ( FractureHandler::FRACTURE_UNCUT * ( M_ID + 1 ) ).add ( M_meshFlat.convex_index() );

	std::ostringstream ss;
	ss << "./matlab/risultati/inizialSaturation_" << M_ID << ".txt";

	std::string name = ss.str();

	std::ofstream exp( name );

	exp << M_ci;

	return;
}// init

void FractureHandler::normalVectorAndMap ( const getfem::mesh_fem& mediumMeshFEMPressure )
{
	// Assegno il vettore normale e la mappa per la frattura
    for ( size_type i = 0; i < M_meshFEM.nb_dof(); ++i )
    {
        const base_node& node = mediumMeshFEMPressure.point_of_basic_dof(i);
        const bgeot::dim_type& dim = M_data.getSpaceDimension();

        scalarVector_Type magnificationMapFactor = M_levelSet->getData()->map_jac( node, dim );

        scalarVector_Type fractureNormal = M_levelSet->getData()->normal_map( node, dim );

        M_magnificationMapFactor1.push_back ( magnificationMapFactor [ 0 ] );
        M_normal1.push_back ( fractureNormal [ 0 ] );
        M_normal2.push_back ( fractureNormal [ 1 ] );
    }

    return;

}// normalVectorAndMap


sizeVector_Type FractureHandler::getDofIntersection( )
{
	sizeVector_Type Dof;

	Dof.push_back( M_data.getInt() );

	return Dof;
}// fetDofIntersection


void FractureHandler::computeInvH ( const BCHandlerPtr_Type& bcHandler )
{
    // Computing h^(-1) on external boudaries of the fracture.
    // Useful for impose the boundary condition with Nitsche penalisation

    const sizeVector_Type& fractureDirichlet = bcHandler->getFractureBC( M_ID )->getDirichlet();
    const size_type shiftFracture = fractureDirichlet.size();

    for ( size_type i = 0; i < shiftFracture; i++ )
    {
        for ( getfem::mr_visitor vis( M_meshFlat.region( fractureDirichlet [ i ] ) ); !vis.finished(); ++vis )
        {
            // Seleziono i grado di libertà corrente
            size_type dof = M_meshFEM.ind_basic_dof_of_element( vis.cv() ) [ 0 ];

            // Stimo la dimensione della mesh
            const scalar_type meshSize = M_meshFlat.convex_radius_estimate( vis.cv() );

            // calcolo h^(-1)
            M_inverseMeshSize [ dof ] = 1.0 / meshSize;

        }
    }

    return;

}// computeInvH


void FractureHandler::setMeshLevelSetFracture ( FractureHandler& otherFracture )
{
    const size_type otherFractureId = otherFracture.getID();
    size_type numIntersect = 0;

    if ( !M_meshLevelSetIntersect[ otherFractureId ].get() )
    {
		M_meshLevelSetIntersect[ otherFractureId ].reset ( new GFMeshLevelSet_Type ( M_meshFlat ) );
        LevelSetHandlerPtr_Type otherLevelSet = otherFracture.getLevelSet();
        M_levelSetIntersect [ otherFractureId ].reset ( new GFLevelSet_Type ( M_meshFlat, 1, false )  );
        M_levelSetIntersect [ otherFractureId ]->reinit();

        const size_type nbDof = M_levelSetIntersect [ otherFractureId ]->get_mesh_fem().nb_basic_dof();

        for ( size_type d = 0; d < nbDof; ++d )
        {
            base_node node = M_levelSetIntersect [ otherFractureId ]->get_mesh_fem().point_of_basic_dof(d);
            base_node mappedNode ( node.size() +1 );

			scalar_type t = d*1./(M_data.getSpatialDiscretization () );
			base_node P (node.size());

			P [0] = t;

			mappedNode [0] = node [0];

			mappedNode [1] = M_levelSet->getData()->y_map ( P );

			M_levelSetIntersect [ otherFractureId ]->values(0)[d] = otherLevelSet->getData()->levelSetFunction ( mappedNode );

        }

    }

    return;

}// setMeshLevelSetFracture


void FractureHandler::updateEtaTangentialInterpolated ( const scalarVector_Type& s )
{
	size_type dof = M_data.getSpatialDiscretization()-1;

	/*
    // Alloco il vettore per M_etaNormalInterpolated
    gmm::resize ( M_etaNormalInterpolated, M_meshFEM.nb_dof() );
    gmm::clear ( M_etaNormalInterpolated );

    std::cout << "   9  " << std::endl;
    // Alloco il vettore per M_etaTangentialInterpolated
    gmm::resize ( M_etaTangentialInterpolated, M_meshFEM.nb_dof() );
    gmm::clear ( M_etaTangentialInterpolated );

	*/

    // Alloco il vettore per M_etaNormalInterpolated
    gmm::resize ( M_etaNormalInterpolated, dof);
    gmm::clear ( M_etaNormalInterpolated );


    // Alloco il vettore per M_etaTangentialInterpolated
    gmm::resize ( M_etaTangentialInterpolated, dof );
    gmm::clear ( M_etaTangentialInterpolated );

    /*
    // Riempio i vettori M_etaNormalInterpolated, M_etaTangentialInterpolated, M_muNormalInterpolated e M_muTangentialInterpolated della frattura
    for ( size_type i = 0; i < M_meshFEM.nb_dof(); ++i )
    {
        M_etaNormalInterpolated [ i ] = M_data.etaNormalDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaNormal()/( s [ i ] );

        M_etaTangentialInterpolated [ i ] = M_data.etaTangentialDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaTangential() /( s [ i ] );
    }
	*/

    scalar_type lambda;



    /*
    // Riempio i vettori M_etaNormalInterpolated, M_etaTangentialInterpolated, M_muNormalInterpolated e M_muTangentialInterpolated della frattura
    for ( size_type i = 0; i < M_meshFEM.nb_dof(); ++i )
    {
    	lambda = M_data.lambdaW( ( M_meshFEM.point_of_dof(i) ) ) + M_data.lambdaNW( ( M_meshFEM.point_of_dof(i) ) );
        M_etaNormalInterpolated [ i ] = M_data.etaNormalDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaNormal()/lambda;

        M_etaTangentialInterpolated [ i ] = M_data.etaTangentialDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaTangential() /lambda;
    }
	*/

    // Riempio i vettori M_etaNormalInterpolated, M_etaTangentialInterpolated, M_muNormalInterpolated e M_muTangentialInterpolated della frattura
    for ( size_type i = 0; i < dof; ++i )
    {
     	lambda = M_data.lambdaW( ( M_meshFEM.point_of_dof(i) ) ) + M_data.lambdaNW( ( M_meshFEM.point_of_dof(i) ) );

    	M_etaNormalInterpolated [ i ] = M_data.etaNormalDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaNormal()/lambda;

        M_etaTangentialInterpolated [ i ] = M_data.etaTangentialDistribution( M_meshFEM.point_of_dof(i) ) * M_data.getEtaTangential() /lambda;
    }

	return;
}// updateEtaTangentialInterpolated
