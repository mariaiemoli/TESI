
#include "../include/MeshHandler.h"

/**************************************************************************/
/*  MeshHandler.cc														  */
/*  Classe che costruisce e manipola la mesh di supporto                  */
/**************************************************************************/

MeshHandler::MeshHandler ( const GetPot& dataFile, const std::string& sectionDomain ):

			M_meshLevelSet ( M_mesh ),
			M_meshType(dataFile( ( sectionDomain + "meshTyper" ).data(), "GT_PK(2,1)" ) ),
			M_spaceDimension(dataFile( ( sectionDomain + "sapceDimension" ).data(), 2.) ),
			M_fEMTypeVector( dataFile ( ( sectionDomain + "fEMTypeVelocity" ).data(), "FEM_RT0(2)" ) ),
			M_meshFEMVector( M_mesh ),
			M_fEMTypeScalar( dataFile ( ( sectionDomain + "fEMPressure" ).data(), "FEM_PK(2,0)" ) ),
			M_meshFEMScalar( M_mesh ),
			M_meshFEMCoefficients( M_mesh ),
			M_integrationTypeVector( dataFile ( ( sectionDomain + "integrationTypeVelocity" ).data(), "IM_TRIANGLE(6)" ) ),
			M_integrationMethodVector( M_mesh ),
			M_integrationTypeScalar( dataFile ( ( sectionDomain + "integrationTypePressure" ).data(), "IM_TRIANGLE(1)" ) ),
			M_integrationMethodScalar( M_mesh ),
			M_meshExternal( dataFile ( ( sectionDomain + "meshExternal" ).data(), "none" ) ),
			M_meshFolder( dataFile ( ( sectionDomain + "meshFolder" ).data(), "./" ) ),
			M_spatialDiscretization( dataFile ( ( sectionDomain + "spatialDiscretization" ).data(), 10 ) ),
			M_inclination( dataFile ( ( sectionDomain + "spatialInclination" ).data(), 0. ) ),
			M_lengthAbscissa( dataFile ( ( sectionDomain + "lengthAbscissa" ).data(), 1. ) ),
			M_lengthOrdinate( dataFile ( ( sectionDomain + "lengthOrdinate").data(), 1. ) ),
			M_lengthQuota( dataFile ( ( sectionDomain + "lengthQuota" ).data(), 1. ) )
{
}// costruttore


void MeshHandler::setUpMesh ( )
{
	if ( M_meshExternal == "none" )
	{
		//------------------M_mediumMesh di Omega--------------------------------
		sizeVector_Type numberSubdivision( M_spaceDimension );
		std::fill( numberSubdivision.begin(), numberSubdivision.end(), M_spatialDiscretization );

        // Geometric transformation usign primal finite elements type
		M_geometricTransformation = bgeot::geometric_trans_descriptor( M_meshType);

		getfem::regular_unit_mesh( M_mesh, numberSubdivision, M_geometricTransformation );

		bgeot::base_matrix transformMatrix( M_spaceDimension, M_spaceDimension );

		scalarVector_Type length(3, 0);
		length [ 0 ] = M_lengthAbscissa;
		length [ 1 ] = M_lengthOrdinate;
		length [ 2 ] = M_lengthQuota;


        for ( size_type i = 0; i < M_spaceDimension; ++i )
        {
            transformMatrix(i, i) = (i < M_spaceDimension) ? length [ i ] : 1.0;
        }
        if ( M_spaceDimension > 1 )
        {
            transformMatrix(0, 1) = M_inclination * M_lengthOrdinate;
        }

        // riscalo la mesh unitaria M_mediumMesh to [M_mediumLengthAbscissa,M_mediumLengthOrdinate]
        M_mesh.transformation(transformMatrix);

	}

	else
		// se M_meshExternal != "none" importo la mesh già costruita
		getfem::import_mesh((M_meshFolder + M_meshExternal).data(), M_mesh);
	
	return;

}// setUpMesh



// Definizione degli elementi finiti
void MeshHandler::setUpFEM ( )
{
	// Dual variable spaces
    getfem::pfem FETypeVelocity = getfem::fem_descriptor(M_fEMTypeVector);

    // Integration type for the dual variable
    getfem::pintegration_method integrationTypeVector = getfem::int_method_descriptor(M_integrationTypeVector);

    // Integration method for the dual variable
    M_integrationMethodVector.set_integration_method(M_mesh.convex_index(), integrationTypeVector);

    // Finite element space for the dual variable
    M_meshFEMVector.set_qdim(M_spaceDimension);
    M_meshFEMVector.set_finite_element(M_mesh.convex_index(), FETypeVelocity);

    // Finite element type for the primal variable
    getfem::pfem FETypePressure = getfem::fem_descriptor(M_fEMTypeScalar);

    // Integration type for the primal variable
    getfem::pintegration_method integrationTypePressure = getfem::int_method_descriptor(M_integrationTypeScalar);

    // Integration method for the primal variable
    M_integrationMethodScalar.set_integration_method(M_mesh.convex_index(), integrationTypePressure);

    //  Finite element space for the primal variable
    M_meshFEMScalar.set_finite_element(M_mesh.convex_index(), FETypePressure);

    // Coefficient: P0 FEM

    // Finite element space the coefficients, the same as the primal finite element space
    M_meshFEMCoefficients.set_finite_element(M_mesh.convex_index(), FETypePressure);

    return;

}// setUpFEM



void MeshHandler::computeMeshMeasures ( )
{
    // Allocate the vector for the circumcenters
    gmm::resize(M_circumcentersAbscissa, M_meshFEMCoefficients.nb_dof());
    gmm::clear(M_circumcentersAbscissa);

    gmm::resize(M_circumcentersOrdinate, M_meshFEMCoefficients.nb_dof());
    gmm::clear(M_circumcentersOrdinate);

    // Allocate the vector for the edge midpoint
    gmm::resize(M_edgeMidpointAbscissa, M_meshFEMVector.nb_dof());
    gmm::clear(M_edgeMidpointAbscissa);

    gmm::resize(M_edgeMidpointOrdinate, M_meshFEMVector.nb_dof());
    gmm::clear(M_edgeMidpointOrdinate);

    // Allocate the vector for the edge length
    gmm::resize(M_edgeLength, M_meshFEMVector.nb_dof());
    gmm::clear(M_edgeLength);
	
	
    // Computes all the circumcenters
	for ( size_type j=0; j< (M_mesh.convex_index()).size(); j++)
    {	
		// Get the coordinates of the vertices of the current triangle
        bgeot::basic_mesh::ref_mesh_pt_ct coordinates = M_mesh.points_of_convex(j);

        for ( size_type i = 0; i < M_mesh.nb_faces_of_convex(j); ++i )
        {
			
            size_type dof = M_meshFEMVector.ind_basic_dof_of_element(j) [ i ];

            M_edgeMidpointAbscissa [ dof ] = 0.5 * (coordinates [ i ] [ 0 ] + coordinates [ (i + 1) % 3 ] [ 0 ]);

            M_edgeMidpointOrdinate [ dof ] = 0.5 * (coordinates [ i ] [ 1 ] + coordinates [ (i + 1) % 3 ] [ 1 ]);
			
            // Computes the edge length
            M_edgeLength [ dof ] = pointDistance(coordinates [ i ] [ 0 ],
									coordinates [ (i + 1) % 3 ] [ 0 ], coordinates [ i ] [ 1 ],
									coordinates [ (i + 1) % 3 ] [ 1 ]);
        }

        // Computes the circumcenter
        for ( size_type i = 0; i < M_mesh.nb_faces_of_convex(j); ++i )
        {
            size_type dof = M_meshFEMVector.ind_basic_dof_of_element(j) [ i ];
			
			M_circumcentersAbscissa [ j ] += 0.3 * M_edgeMidpointAbscissa [ dof ];
			
            M_circumcentersOrdinate [ j ] += 0.3 * M_edgeMidpointOrdinate [ dof ];
        }
	
    }
	
    gmm::resize(M_circumcentersDistance, M_meshFEMVector.nb_dof());
    gmm::clear(M_circumcentersDistance);

    // Computes all the distances between the circumcenters of adjacent triangels
    for ( size_type j=0; j< (M_mesh.convex_index()).size(); j++)
	{
		for ( size_type i = 0; i < M_mesh.nb_faces_of_convex(j); ++i )
        {
            // For each face get the corresponding neighbour
			size_type neighbour = M_mesh.neighbour_of_convex(j, i);

            size_type dof = M_meshFEMVector.ind_basic_dof_of_element(j) [ i ];
 
            // Check if the neighbour exist, i.e. not a boundary element
            if ( neighbour != size_type(-1) )
            {

                // Check if the position in M_circumcentersDistance is jet filled
                if ( M_circumcentersDistance [ dof ] == 0 )
                {
                    // Computes the distance of the circumcenters
                    M_circumcentersDistance [ dof ] = pointDistance( M_circumcentersAbscissa [ neighbour ],
                    												 M_circumcentersAbscissa [ j ],
																	 M_circumcentersOrdinate [ neighbour ],
																	 M_circumcentersOrdinate [ j ]);

                }
            }
            else
            {
                // Boundary face, the position is not jet filled
                M_circumcentersDistance [ dof ] = pointDistance( M_edgeMidpointAbscissa [ dof ],
																 M_circumcentersAbscissa [ j ],
																 M_edgeMidpointOrdinate [ dof ],
																 M_circumcentersOrdinate [ j ]);
            }
        }
    }

    // Computing h^(-1) on external boudaries.
    // Useful for impose the boundary condition with Nitsche penalisation

    // Allocate the vector for the M_mediumMesh size
    gmm::resize(M_meshSize, M_meshFEMScalar.nb_dof());
    gmm::clear(M_meshSize);

    // Allocate the vector for the inverse of the M_mediumMesh size
    gmm::resize(M_inverseMeshSize, M_meshFEMScalar.nb_dof());
    gmm::clear(M_inverseMeshSize);

    for ( size_type j=0; j< (M_mesh.convex_index()).size(); j++)
    {
        // Select the current dof
        size_type dofScalar = M_meshFEMScalar.ind_basic_dof_of_element(j) [ 0 ];

        scalarVector_Type edges(3, 0.);

		for ( size_type i = 0; i < M_mesh.nb_faces_of_convex(j); ++i )
        {

            size_type dofVector = M_meshFEMVector.ind_basic_dof_of_element( j ) [ i ];

            // Computes the edge length
            edges [ i ] = M_edgeLength [ dofVector ];

        }

        // Estimate the element size
        M_meshSize [ dofScalar ] = *(std::max_element ( edges.begin(), edges.end() ));

        // Compute h^(-1)
        M_inverseMeshSize [ dofScalar ] = 1.0 / M_meshSize [ dofScalar ];
    }

    return;
    
}// computeMeshMeasures
