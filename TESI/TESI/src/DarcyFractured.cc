
#include "../include/DarcyFractured.h"

/**************************************************************************/
/*  DarcyFractured.cc													  */
/*  Classe che assembla e risolve il sistema per il problema di Darcy	  */ 
/*  sulle fratture											              */
/**************************************************************************/


DarcyFractured::DarcyFractured ( const MediumDataPtr_Type& medium,
                                 const MeshHandlerPtr_Type& mesh,
                                 const BCHandlerPtr_Type& bcHandler,
                                 const FracturesSetPtr_Type& fractures,
                                 const ExporterPtr_Type& exporter ) :
								M_mediumData ( medium ),
								M_mesh ( mesh ),
								M_bcHandler ( bcHandler ),
								M_fractures ( fractures ),
								M_exporter ( exporter ),
								M_fractureEtaNormalOnMedium ( M_fractures->getNumberFractures () ),
								M_fractureVelocity ( M_fractures->getNumberFractures () ),
								M_fracturePressure ( M_fractures->getNumberFractures () )
{} // costruttore


void DarcyFractured::init ( )
{

    // Alloco i dati
    const size_type shifCoefficients = M_mesh->getMeshFEMCoefficients().nb_dof();

    const size_type numberFractures = M_fractures->getNumberFractures();

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Alloco il vettore per M_fracturesEtaNormalOnMedium
        gmm::resize(M_fractureEtaNormalOnMedium [ f ], shifCoefficients);
        gmm::clear(M_fractureEtaNormalOnMedium [ f ]);

    }

    return;

} // init


void DarcyFractured::assembly ( const GetPot& dataFile )
{
	// numero totale delle fratture
    const scalar_type numberFractures = M_fractures->getNumberFractures();

    // vettore dei gradi di libertà per la pressione per le fratture
    sizeVector_Type fractureNumberDOFPressure(numberFractures);

    // numero di gradi di libertà totali per la pressione (con quelli estesi) per ogni frattura
    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    // vettore dei gradi di libertà per la velocità per le fratture
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);

    // numero di gradi di libertà totali per la velocità (con quelli estesi) per ogni frattura
    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    // numero totale dei gradi di libertà ( velocità + pressione ) per ogni frattura
    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    // numero di gradi di libertà vincolati alla condizione al bordo per ogni frattura
    sizeVector_Type fractureNumberBoundaryDOF(numberFractures);

    // numero complessivo dei gradi di libertà ( pressione + velocità )
    size_type fractureTotalNumberDOFVelocityPressure = 0;

    // numero complessivo dei gradi di libertà di pressione
    size_type fractureTotalNumberDOFPressure = 0;

    // numero totale di intersezioni
    size_type fractureNumberBifurcation = 0;
    size_type fractureNumberCross = 0;
    size_type globalFractureNumber = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	// numero di gradi di libertà per la pressione senza quelli estesi per la frattura f
        fractureNumberDOFPressure [ f ] = M_fractures->getFracture( f )->getMeshFEM().nb_dof();

    	// numero di gradi di libertà per la pressione con quelli estesi per la frattura f
        fractureNumberGlobalDOFPressure [ f ] = fractureNumberDOFPressure [ f ];

    	// numero totale di gradi di libertà per la pressione
        fractureTotalNumberDOFPressure += fractureNumberGlobalDOFPressure [ f ];

    	// numero di gradi di libertà per la velocità senza quelli estesi per la frattura f
        fractureNumberDOFVelocity [ f ] = M_fractures->getFracture( f )->getMeshFEM2().nb_dof();

    	// numero di gradi di libertà per la velocità con quelli estesi per la frattura f
        fractureNumberGlobalDOFVelocity [ f ] = fractureNumberDOFVelocity [ f ];

    	// numero totale di gradi di libertà per la velocità e la pressione per la frattura f
        fractureNumberDOFVelocityPressure [ f ] = fractureNumberGlobalDOFVelocity [ f ] + fractureNumberGlobalDOFPressure [ f ];

    	// numero totale di gradi di libertà per la pressione e la velocità
        fractureTotalNumberDOFVelocityPressure += fractureNumberDOFVelocityPressure [ f ];

    	// numero di gradi di libertà per i bordi per la frattura f
        fractureNumberBoundaryDOF [ f ] = M_bcHandler->getFractureBC(f)->getMeshFEM().nb_dof();

    }

    // Numero intersezioni

    fractureNumberBifurcation = M_fractures->getIntersections ()->getNumberBifurcation ();
    fractureNumberCross = M_fractures->getIntersections ()->getNumberCross ();

    globalFractureNumber = fractureNumberBifurcation + fractureNumberCross;


    // Inizializziamo tutte le matrici a blocchi, la matrice globale e il termine di destra per il sistema e il vettore delle soluzioni

    // Allochiamo la matrice globale del sistema:  M_darcyGlobalMatrix
    M_globalMatrix.reset( new sparseMatrix_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber,
             	 	 	 	 	 	 	 	 	 fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));

    gmm::clear( *M_globalMatrix );

    // Allochiamo il vettore del termine noto di destra del sistema: M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset( new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear( *M_globalRightHandSide );

    // Allochiamo il vettore globale del termine incognito del sistema: M_darcyVelocityAndPressure
    M_velocityAndPressure.reset(new scalarVector_Type( fractureTotalNumberDOFVelocityPressure + globalFractureNumber ));
    gmm::clear(*M_velocityAndPressure);


    // Matrici a blocchi per le fratture
    sparseMatrixPtrContainer_Type A11F(numberFractures), A12F(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Alloco la matrice A11F per la fratture f-esima
        A11F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFVelocity [ f ]) );
        gmm::clear(*(A11F [ f ]));

        // Alloco la matrice A12F per la fratture f-esima
        A12F [ f ].reset(new sparseMatrix_Type ( fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ]) );
        gmm::clear(*(A12F [ f ]));
    }

    sizeVector_Type shiftIntersect ( numberFractures );

    shiftIntersect [ 0 ] = 0;

    // Calcolo i blocchi di matrici per ogni frattura
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        std::cout <<std::endl<< "Fracture " << f << std::endl;

        // Computes the matrix \int_\gamma \tau_i \cdot \tau_j
        getfem::darcy_A11F ( A11F [ f ], M_fractures->getFracture( f ),
                             M_mediumData->getPenaltyVector(),
                             M_fractures->getFracture( f )->getEtaTangentialInterpolated(),
                             M_bcHandler->getFractureBC(f)->getDirichlet(), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        // Computes the matrix \int_\gamma \nabla v_i \cdot \tau_j
        getfem::darcy_A12F ( A12F [ f ], M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

        if( f != 0)
        {
        	shiftIntersect [ f ] = shiftIntersect [ f-1 ] + fractureNumberDOFVelocityPressure [ f-1 ];
        }

    }

    std::cout << std::endl;


    // Introduco il termine legato all'intersezione
    // A seconda del tipo di intersezione, Cross o Biforcazione, devo introdurre dei termini diversi

    IntersectDataContainer_Type IntBifurcation = M_fractures->getIntersections ()-> getBifurcationIntersections ();

    IntersectDataContainer_Type IntCross = M_fractures->getIntersections ()-> getCrossIntersections ();

        ///////////////////////////////////////////////////////////////////////

	/*
	* 		Copy blocks into the system matrix M_darcyGlobalMatrix
	*
	* 				[  A11  A12  0          ]
	* 		    M = [ -A12  A22  0          ]
	* 		    	[  0    0    A11F  A12F ]
	* 		        [  0    0   -A12F  0    ]
	*/

        ///////////////////////////////////////////////////////////////////////

	// Shift for the fracture
	size_type fractureShift = 0;

	for ( size_type f = 0; f < numberFractures; ++f )
	{
	   // Copy the matrix A11F in M_darcyGlobalMatrix in the correct position
	   gmm::copy( *( A11F [ f ] ), gmm::sub_matrix( *M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ] ),
	                   										     gmm::sub_interval(fractureShift, fractureNumberGlobalDOFVelocity [ f ] )));

	   // Copy the matrix A12F in M_darcyGlobalMatrix in the correct position
	   gmm::copy( *( A12F [ f ] ), gmm::sub_matrix(*M_globalMatrix, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ]),
	   														  gmm::sub_interval(fractureShift + fractureNumberGlobalDOFVelocity [ f ],
	   																  	        fractureNumberGlobalDOFPressure [ f ])));

	   // Copy the matrix -A12F in M_darcyGlobalMatrix in the correct position
	   gmm::copy(gmm::transposed(gmm::scaled(*(A12F [ f ]), -1.0)),
	   						  gmm::sub_matrix( *M_globalMatrix,
	   								  gmm::sub_interval( fractureShift + fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ]),
	   								  gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ])));

	   // Update the shift
	   fractureShift += fractureNumberDOFVelocityPressure [ f ];

	}

	sizeVector_Type shiftVelocity ( numberFractures );
	shiftVelocity [ 0 ] = fractureNumberDOFVelocity [ 0 ];

	for ( size_type f = 0; f < numberFractures; ++f )
	{
		if( f != 0)
		{
			shiftVelocity [ f ] = shiftVelocity [ f-1 ] + fractureNumberGlobalDOFPressure [ f-1 ] + fractureNumberGlobalDOFVelocity [ f ];
		}

	}

    // Aggiorno ora la matrice globale imponendo le condizioni di interfaccia per la biforcazione

    for ( size_type i = 0; i < IntBifurcation.size(); i++ )
    {
    	sparseMatrixPtr_Type Aup0, Aup1, Aup2, Aup3;

    	FractureHandlerPtr_Type f0 = IntBifurcation [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntBifurcation [ i ].getFracture (1);
    	FractureHandlerPtr_Type f2 = IntBifurcation [ i ].getFracture (2);

    	size_type id0 = f0->getID();
    	size_type id1 = f1->getID();
    	size_type id2 = f2->getID();

        std::cout << "Bifurcation: " << std::endl;
		std::cout << " Fractures " << id0 << ", " << id1 << ", " << id2 << std::endl;

		//const size_type globalIndex = GlobalIndex_Bifurcation( M_fractures , id0,id1, id2 );
		const size_type globalIndex = i;
		const size_type Index = fractureTotalNumberDOFVelocityPressure + globalIndex;

        Aup0.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup0);

        Aup1.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup1);

        Aup2.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup2);

        Aup3.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup3);


        MatrixBifurcationHandler_Type Matrix( dataFile );
		FracturePtrContainer_Type Fracture( 3 );
		Fracture[ 0 ] =f0;
		Fracture[ 1 ] =f1;
		Fracture[ 2 ] =f2;


		FracturePtrContainer_Type Fracture_copy( 3 );
		Fracture_copy = Fracture;

		Matrix.setMatrices( Fracture_copy );

		Matrix4d T = Matrix.T();

		scalar_type s = 0.;
		Matrix.computeScap ( s );

		scalarVector_Type DOF( 3 );
		scalarVector_Type DOF_v( 3 );

		for ( size_type j = 0; j < 3; j++ )
		{
			DOF [ j ] = Fracture [ j ] -> getData().getInt();

			if( DOF[ j ] == 0 )
			{
				DOF_v[ j ] = DOF[ j ];
			}
			else
			{
				DOF_v[ j ] = DOF[ j ] + 1.;
			}
		}


		getfem::setAup_i( Aup0, 0 , id0, id1, id2, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup1, 1 , id1, id0, id2, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup2, 2 , id2, id0, id1, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup3, 3 , id2, id0, id1, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index, s );

		gmm::copy(*Aup0, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id0 ] + DOF[ 0 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup1, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id1 ] + DOF[ 1 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup2, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id2 ] + DOF[ 2 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup3, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure + fractureNumberCross + i, 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));


    }

    // Aggiorno ora la matrice globale imponendo le condizioni di interfaccia per il cross

    for ( size_type i = 0; i < IntCross.size(); i++ )
    {
    	sparseMatrixPtr_Type Aup0, Aup1, Aup2, Aup3, Aup4;

    	FractureHandlerPtr_Type f0 = IntCross [ i ].getFracture (0);
    	FractureHandlerPtr_Type f1 = IntCross [ i ].getFracture (1);
    	FractureHandlerPtr_Type f2 = IntCross [ i ].getFracture (2);
    	FractureHandlerPtr_Type f3 = IntCross [ i ].getFracture (3);

    	size_type id0 = f0->getID();
    	size_type id1 = f1->getID();
    	size_type id2 = f2->getID();
    	size_type id3 = f3->getID();

        std::cout << "Cross: " << std::endl;
		std::cout << " Fractures " << id0 << ", " << id1 << ", " << id2 << ", " << id3 << std::endl;


		const size_type globalIndex = i;
		const size_type Index = fractureTotalNumberDOFVelocityPressure + globalIndex;

        Aup0.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup0);

        Aup1.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup1);

        Aup2.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
        gmm::clear(*Aup2);

        Aup3.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup3);

        Aup4.reset ( new sparseMatrix_Type ( 1, fractureTotalNumberDOFVelocityPressure + globalFractureNumber) );
	    gmm::clear(*Aup4);


        MatrixBifurcationHandler_Type Matrix( dataFile, "Cross" );
		FracturePtrContainer_Type Fracture( 4 );
		Fracture[ 0 ] =f0;
		Fracture[ 1 ] =f1;
		Fracture[ 2 ] =f2;
		Fracture[ 3 ] =f3;


		FracturePtrContainer_Type Fracture_copy( 4 );
		Fracture_copy = Fracture;

		Matrix.setMatrices( Fracture_copy );

		Matrix4d T = Matrix.T();

		scalar_type s = 0.;
		Matrix.computeScap ( s );

		scalarVector_Type DOF( 4 );
		scalarVector_Type DOF_v( 4 );

		sizeVector_Type ID ( 4 );


		for ( size_type j = 0; j < 4; j++ )
		{
			DOF [ j ] = Fracture [ j ] -> getData().getInt();

			if( DOF[ j ] == 0 )
			{
				DOF_v[ j ] = DOF[ j ];
			}
			else
			{
				DOF_v[ j ] = DOF[ j ] + 1.;
			}

			ID [ j ] = IntCross [ i ].getFracture (j)->getID();
		}


		getfem::setAup_i( Aup0, 0 , ID, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup1, 1 , ID, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup2, 2 , ID, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup2, 3 , ID, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index );
		getfem::setAup_i( Aup3, 4 , ID, DOF, DOF_v, shiftIntersect, fractureNumberGlobalDOFVelocity,  T, Index, s );


		gmm::copy(*Aup0, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id0 ] + DOF[ 0 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup1, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id1 ] + DOF[ 1 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup2, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id2 ] + DOF[ 2 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

		gmm::copy(*Aup3, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( shiftVelocity[ id3 ] + DOF[ 3 ], 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));


		gmm::copy(*Aup4, gmm::sub_matrix(*M_globalMatrix,
	    		gmm::sub_interval( fractureTotalNumberDOFVelocityPressure + i, 1),
	    		gmm::sub_interval( 0, fractureTotalNumberDOFVelocityPressure + globalFractureNumber ) ));

    }

    //Costruiamo il termine noto

    // Vector, Neumann conditions
    scalarVectorPtrContainer_Type BstressF(numberFractures);

    // Pressure term (Neumann)
    scalarVectorPtrContainer_Type PneumannF(numberFractures);

    // Velocity term (Dirichilet)
    scalarVectorPtrContainer_Type VdirichletF(numberFractures);

    // Right Hand Side: velocity and pressure
    scalarVectorPtrContainer_Type B_vF(numberFractures);

    // Right Hand Side: velocity and pressure
    scalarVectorPtrContainer_Type B_pF(numberFractures);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        BstressF [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFVelocity [ f ], 0.));

        PneumannF [ f ].reset(new scalarVector_Type( fractureNumberBoundaryDOF [ f ], 0.));

        VdirichletF [ f ].reset(new scalarVector_Type( fractureNumberDOFPressure [ f ], 0.));

        B_vF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFVelocity [ f ], 0.));

        B_pF [ f ].reset(new scalarVector_Type(fractureNumberGlobalDOFPressure [ f ], 0.));

    }

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        const BCPtr_Type& fractureBC = M_bcHandler->getFractureBC( f );

        for ( size_type i = 0; i < fractureNumberBoundaryDOF [ f ]; i++ )
        {
            const base_node node = fractureBC->getMeshFEM().point_of_basic_dof( i );

            (*(PneumannF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().pressureExact( node );
        }

        (*(PneumannF [ f ])) [ 0 ] *= -1;
    }

    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	std::cout << "Fracture " << f << std::endl;

        // Computes the right hand side rhs which stores the boundary conditions for the fracture
        getfem::darcy_dataF ( BstressF [ f ], B_vF [ f ], M_bcHandler,
                              M_fractures->getFracture( f ), M_mediumData->getPenaltyVector(),
                              M_mediumData->getInvK(), PneumannF [ f ], VdirichletF [ f ] );

        std::cout<< std::endl;
    }

    fractureShift = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        for ( size_type i = 0; i < fractureNumberGlobalDOFVelocity [ f ]; ++i )
        {
            (*M_globalRightHandSide) [ fractureShift + i ] += (*(BstressF [ f ])) [ i ];

            (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_vF [ f ])) [ i ];
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }


    scalarVectorPtrContainer_Type divF(numberFractures);

    fractureShift = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        divF [ f ].reset(new scalarVector_Type(fractureNumberDOFPressure [ f ], 0.));

        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; i++ )
        {
            const base_node node = M_fractures->getFracture( f )->getMeshFEM().point_of_basic_dof( i );

            (*(divF [ f ])) [ i ] = M_fractures->getFracture( f )->getData().darcySource( node );
        }

        gmm::clear(*(B_pF [ f ]));

        getfem::assembling_Source_BoundaryF ( B_pF [ f ], divF [ f ], M_fractures->getFracture( f ), FractureHandler::FRACTURE_UNCUT * ( f + 1 ) );

    }

    fractureShift = 0;

    for ( size_type f = 0; f < numberFractures; ++f )
    {

        for ( size_type i = fractureNumberGlobalDOFVelocity [ f ]; i < fractureNumberDOFVelocityPressure [ f ]; ++i )
        {
                (*M_globalRightHandSide) [ fractureShift + i ] += (*(B_pF [ f ])) [ i - fractureNumberGlobalDOFVelocity [ f ] ];
        }

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];

    }

    // Curo la matrice
    for ( size_type i = 0; i < M_globalRightHandSide->size(); ++i )
    {
        scalar_type somma = 0;

        for ( size_type j = 0; j < M_globalRightHandSide->size(); ++j )
        {
            somma += std::fabs( (*M_globalMatrix) (i, j) );
        }

        if ( somma == 0 )
        {
            (*M_globalRightHandSide) [ i ] = 0.;
            (*M_globalMatrix) (i, i) = 1.;
        }
    }

    M_exporter->spy(M_globalMatrix, "./Matlab/matrice.mm");
    M_exporter->spy(M_globalRightHandSide, "./Matlab/rhs.mm");

    return;

}// assembly


// Solve the Darcy for the governing flux and do the time loop for the evolution problem
void DarcyFractured::solve ( )
{
	// numero totale delle fratture
	const scalar_type numberFractures = M_fractures->getNumberFractures();

	// vettore dei gradi di libertà per la velocità per le fratture
    sizeVector_Type fractureNumberDOFVelocity(numberFractures);

    // numero di gradi di libertà totali per la velocità (con quelli estesi) per ogni frattura
    sizeVector_Type fractureNumberGlobalDOFVelocity(numberFractures);

    // vettore dei gradi di libertà per la pressione per le fratture
    sizeVector_Type fractureNumberDOFPressure(numberFractures);

    // numero di gradi di libertà totali per la pressione (con quelli estesi) per ogni frattura
    sizeVector_Type fractureNumberGlobalDOFPressure(numberFractures);

    // numero totale dei gradi di libertà ( velocità + pressione ) per ogni frattura
    sizeVector_Type fractureNumberDOFVelocityPressure(numberFractures);

    // numero complessivo dei gradi di libertà ( pressione + velocità )
    size_type fractureTotalNumberDOFVelocityPressure(0);


    for ( size_type f = 0; f < numberFractures; ++f )
    {
    	// numero di gradi di libertà per la velocità senza quelli estesi per la frattura f
        fractureNumberDOFVelocity [ f ] = M_fractures->getFracture( f )->getMeshFEM2().nb_dof();

        // numero di gradi di libertà per la velocità con quelli estesi per la frattura f
        fractureNumberGlobalDOFVelocity [ f ] = fractureNumberDOFVelocity [ f ];

        // numero di gradi di libertà per la pressione senza quelli estesi per la frattura f
        fractureNumberDOFPressure [ f ] = M_fractures->getFracture( f )->getMeshFEM().nb_dof();

        // numero di gradi di libertà per la pressione con quelli estesi per la frattura f
        fractureNumberGlobalDOFPressure [ f ] = fractureNumberDOFPressure [ f ];

        // numero totale di gradi di libertà per la velocità e la pressione per la frattura f
        fractureNumberDOFVelocityPressure [ f ] = fractureNumberGlobalDOFVelocity [ f ] + fractureNumberGlobalDOFPressure [ f ];

        // numero totale di gradi di libertà per la pressione e la velocità
        fractureTotalNumberDOFVelocityPressure += fractureNumberDOFVelocityPressure [ f ];

        M_fracturePressure [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFPressure [ f ], 0));

        M_fractureVelocity [ f ].reset(new scalarVector_Type( fractureNumberGlobalDOFVelocity [ f ], 0));

    }

    // Solve the Darcy problem
    std::cout << std::endl;
    std::cout << std::endl << "Solving problem in Omega..." << std::flush;
    scalar_type roundConditionNumber;

    SuperLU_solve(*M_globalMatrix, *M_velocityAndPressure,
                  *M_globalRightHandSide, roundConditionNumber);

    size_type fractureShift = 0;

    size_type NumBifurcation = M_fractures->getIntersections ()->getNumberBifurcation ();

    size_type NumIntersection = NumBifurcation;

    scalarVectorPtr_Type M_intersectionPressure;
    M_intersectionPressure.reset(new scalarVector_Type( NumIntersection, 0));

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        // Extract the dual in the fracture
        gmm::copy( gmm::sub_vector( *M_velocityAndPressure, gmm::sub_interval( fractureShift, fractureNumberGlobalDOFVelocity [ f ])),
                  *( M_fractureVelocity [ f ] ) );

        // Extract the primal in the fracture
        gmm::copy( gmm::sub_vector( *M_velocityAndPressure,
        							gmm::sub_interval( fractureShift + fractureNumberGlobalDOFVelocity [ f ], fractureNumberGlobalDOFPressure [ f ])),
        		   *(M_fracturePressure [ f ]));

        // Update the shift
        fractureShift += fractureNumberDOFVelocityPressure [ f ];
    }

    // Extract the pressure at the intersections
    gmm::copy( gmm::sub_vector( *M_velocityAndPressure,
    							gmm::sub_interval( fractureShift, NumIntersection )),
    		   *( M_intersectionPressure ));


    // the extra dof start at the end of the fracture dofs

    getfem::pfem fractureFETypePressure;
    getfem::pfem fractureFETypeVelocity;

    if ( numberFractures > 0 )
    {
        fractureFETypePressure = getfem::fem_descriptor( M_fractures->getFracture( 0 )->getData().getFEMType());
        fractureFETypeVelocity = getfem::fem_descriptor( M_fractures->getFracture( 0 )->getData().getFEMType2());
    }


    // 		********************************************		//
    //--------export solution in the fracture-------------------
    // 		********************************************		//


    for ( size_type f = 0; f < numberFractures; ++f )
    {
        const FractureHandlerPtr_Type& fracture = M_fractures->getFracture(f);
        std::ostringstream osFileName;


        getfem::mesh meshFcutFlat = M_fractures->getFracture( f )-> getMeshFlat ();
        getfem::mesh meshFcutMapped;


        getfem::mesh_fem meshFEMcutFlat ( meshFcutFlat, fracture->getMeshFEM().get_qdim());
        meshFEMcutFlat.set_finite_element ( fracture->getMeshFEM2().fem_of_element(0) );

        scalarVector_Type ordinataUncut ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type ordinataCut ( meshFEMcutFlat.nb_dof(), 0. );

        for ( size_type i = 0; i < fractureNumberDOFVelocity [ f ]; ++i )
        {
            const bgeot::base_node P = fracture->getMeshFEM2().point_of_basic_dof(i);

            bgeot::base_node P1 = fracture->getMesh().points( ) [i];


            scalar_type c= 1./fracture->getData().getSpatialDiscretization ();

            P1 [0] =  i*c;

            ordinataUncut [ i ] = fracture->getLevelSet()->getData()->y_map(P1);
        }

        getfem::interpolation ( fracture->getMeshFEM2(), meshFEMcutFlat, ordinataUncut, ordinataCut );

        for ( size_type i = 0; i < meshFEMcutFlat.nb_dof(); ++i )
        {
            bgeot::base_node P ( fracture->getData().getSpaceDimension() + 2 );

            P [ 0 ] = meshFEMcutFlat.point_of_basic_dof(i)[0];
            P [ 1 ] = ordinataCut [ i ];

            meshFcutMapped.add_point(P);
        }

        const size_type numConvex = meshFcutFlat.convex_index().size();

        for ( size_type i = 0; i < numConvex; ++i )
        {
            std::vector<bgeot::size_type> point(3);

            point [ 0 ] = meshFcutFlat.ind_points_of_convex(i) [ 0 ];
            point [ 1 ] = meshFcutFlat.ind_points_of_convex(i) [ 1 ];
            point [ 2 ] = meshFcutFlat.ind_points_of_convex(i) [ 2 ];

            meshFcutMapped.add_convex( fracture->getGeometricTransformation(), point.begin() );
        }



        // per ogni frattura esporto la mesh reale
        osFileName << "cmeshF2-" << f << ".vtk";
        exportMesh(M_exporter->getFolder() + osFileName.str().c_str(), meshFcutMapped );

        getfem::mesh_fem mfproj ( meshFcutMapped, fracture->getMeshFEM().get_qdim() );
        getfem::mesh_fem mfproj_v ( meshFcutMapped, fracture->getMeshFEM2().get_qdim() );

        mfproj.set_classical_discontinuous_finite_element ( 0, 0.01 );
        mfproj_v.set_classical_discontinuous_finite_element ( 0, 0.01 );

        scalarVector_Type fracturePressureMeanUNCUT ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureIntersection ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureInCut ( fractureNumberDOFPressure [ f ], 0. );
        scalarVector_Type fracturePressureOutCut ( fractureNumberDOFPressure [ f ], 0. );

        scalarVector_Type fractureVelocityMeanUNCUT ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type fractureVelocityInCut ( fractureNumberDOFVelocity [ f ], 0. );
        scalarVector_Type fractureVelocityOutCut ( fractureNumberDOFVelocity [ f ], 0. );


        for ( size_type i = 0; i < fractureNumberDOFPressure [ f ]; ++i )
        {
            const size_type el =  fracture->getMeshFEM().first_convex_of_basic_dof(i);

            if ( fracture->getMeshFlat().region ( FractureHandler::FRACTURE_UNCUT * ( f + 1 ) ).is_in(el) )
            {
                fracturePressureMeanUNCUT [ i ] = (*(M_fracturePressure [ f ])) [ i ];
            }
            else	//	************************	//
            {
            	fracturePressureIntersection [ i ] = (*(M_fracturePressure [ f ])) [ i ];
            }
        }

        for ( size_type i = 0; i < fractureNumberDOFVelocity [ f ]; ++i )
        {
            const size_type el =  fracture->getMeshFEM2().first_convex_of_basic_dof(i);

            if ( fracture->getMeshFlat().region ( FractureHandler::FRACTURE_UNCUT * ( f + 1 ) ).is_in(el) )
            {
                fractureVelocityMeanUNCUT [ i ] = (*(M_fractureVelocity [ f ])) [ i ];
            }
        }

        /*
         * getfem risolve il sistema per la mesh " piatta ", per avere i corretti valori di pressione devo interpolare
         * sulla mesh mappata i valori che ho ottenuto
         *
         */
        scalarVector_Type fracturePressureMeanUNCUTInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanIntersectionInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanInCutInterpolated ( mfproj.nb_dof(), 0. );
        scalarVector_Type fracturePressureMeanOutCutInterpolated ( mfproj.nb_dof(), 0. );

        scalarVector_Type fractureVelocityMeanUNCUTInterpolated ( mfproj_v.nb_dof(), 0. );
        scalarVector_Type fractureVelocityMeanInCutInterpolated ( mfproj_v.nb_dof(), 0. );
        scalarVector_Type fractureVelocityMeanOutCutInterpolated ( mfproj_v.nb_dof(), 0. );


        getfem::mesh_fem mfprojUncut ( fracture->getMesh(), fracture->getMeshFEM().get_qdim());
        mfprojUncut.set_finite_element ( fractureFETypePressure );

        getfem::mesh_fem mfprojUncut_v ( fracture->getMesh(), fracture->getMeshFEM2().get_qdim());
        mfprojUncut_v.set_finite_element ( fractureFETypeVelocity );

        getfem::interpolation ( mfprojUncut, mfproj, fracturePressureMeanUNCUT, fracturePressureMeanUNCUTInterpolated );


        getfem::interpolation ( mfprojUncut, mfproj, fracturePressureIntersection, fracturePressureMeanIntersectionInterpolated );


        getfem::interpolation ( mfprojUncut_v, mfproj_v, fractureVelocityMeanUNCUT, fractureVelocityMeanUNCUTInterpolated );


        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( fracture->getLevelSetIntersect( otherFracture ).get() )
            {
                getfem::mesh_region regionMesh = fracture->getMeshFlat().region ( FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
                size_type i_cv = 0;
                dal::bit_vector bc_cv = regionMesh.index();
                gmm::clear ( fracturePressureInCut );
                gmm::clear ( fracturePressureOutCut );

                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    const size_type ibase = fracture->getMeshFEM().ind_basic_dof_of_element ( i_cv )[0];

                    const size_type position = 0.;

                    const base_node point = mfprojUncut.point_of_basic_dof ( ibase );
                    const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->ylevelSetFunction( point );

                    if ( levelSetValue < 0 )
                    {
                        fracturePressureInCut [ ibase ] = (*(M_fracturePressure [ f ])) [ ibase ];
                        fracturePressureOutCut [ ibase ] = (*(M_fracturePressure [ f ])) [ fractureNumberDOFPressure [ f ] + position ];
                    }
                    else
                    {
                        fracturePressureOutCut [ ibase ] = (*(M_fracturePressure [ f ])) [ ibase ];
                        fracturePressureInCut [ ibase ] = (*(M_fracturePressure [ f ])) [ fractureNumberDOFPressure [ f ] + position ];
                    }

                }

                gmm::clear ( fracturePressureMeanInCutInterpolated );

                getfem::interpolation ( mfprojUncut, mfproj,
                                        fracturePressureInCut, fracturePressureMeanInCutInterpolated );

                gmm::clear ( fracturePressureMeanOutCutInterpolated );

                getfem::interpolation ( mfprojUncut, mfproj,
                                        fracturePressureOutCut, fracturePressureMeanOutCutInterpolated );

                i_cv = 0;
                bc_cv = meshFcutMapped.convex_index();
                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    getfem::mesh_fem::ind_dof_ct idofs = mfproj.ind_basic_dof_of_element ( i_cv );
                    for ( size_type k = 0; k < idofs.size(); ++k )
                    {
                        const base_node node = mfproj.point_of_basic_dof(idofs [ k ]);
                        const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->ylevelSetFunction( node );
                        if ( levelSetValue < 0 )
                        {
                            fracturePressureMeanUNCUTInterpolated [ idofs[k] ] += fracturePressureMeanInCutInterpolated [ idofs[k] ];
                        }
                        else
                        {
                            fracturePressureMeanUNCUTInterpolated [ idofs[k] ] += fracturePressureMeanOutCutInterpolated [ idofs[k] ];
                        }

                    }

                }

            }
        }

        for ( size_type otherFracture = 0; otherFracture < numberFractures; ++otherFracture )
        {
            if ( fracture->getLevelSetIntersect( otherFracture ).get() )
            {
                getfem::mesh_region regionMesh = fracture->getMeshFlat().region ( FractureHandler::FRACTURE_INTERSECT * ( f + 1 ) + otherFracture + 1 );
                size_type i_cv = 0;
                dal::bit_vector bc_cv = regionMesh.index();
                gmm::clear ( fractureVelocityInCut );
                gmm::clear ( fractureVelocityOutCut );

                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    const size_type ibase = fracture->getMeshFEM2().ind_basic_dof_of_element ( i_cv )[0];

                    const size_type position = 0.;

                    const base_node point = mfprojUncut_v.point_of_basic_dof ( ibase );
                    const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->ylevelSetFunction( point );

                    if ( levelSetValue < 0 )
                    {
                        fractureVelocityInCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ ibase ];
                        fractureVelocityOutCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ fractureNumberDOFVelocity [ f ] + position ];
                    }
                    else
                    {
                        fractureVelocityOutCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ ibase ];
                        fractureVelocityInCut [ ibase ] = (*(M_fractureVelocity [ f ])) [ fractureNumberDOFVelocity [ f ] + position ];
                    }

                }

                gmm::clear ( fractureVelocityMeanInCutInterpolated );

                getfem::interpolation ( mfprojUncut_v, mfproj_v,
                                        fractureVelocityInCut, fractureVelocityMeanInCutInterpolated );

                gmm::clear ( fractureVelocityMeanOutCutInterpolated );

                getfem::interpolation ( mfprojUncut_v, mfproj_v,
                                        fractureVelocityOutCut, fractureVelocityMeanOutCutInterpolated );

                i_cv = 0;
                bc_cv = meshFcutMapped.convex_index();
                for ( i_cv << bc_cv; i_cv != size_type(-1); i_cv << bc_cv )
                {
                    getfem::mesh_fem::ind_dof_ct idofs = mfproj_v.ind_basic_dof_of_element ( i_cv );
                    for ( size_type k = 0; k < idofs.size(); ++k )
                    {
                        const base_node node = mfproj_v.point_of_basic_dof(idofs [ k ]);
                        const scalar_type levelSetValue = M_fractures->getFracture( otherFracture )->getLevelSet()->getData()->ylevelSetFunction( node );
                        if ( levelSetValue < 0 )
                        {
                            fractureVelocityMeanUNCUTInterpolated [ idofs[k] ] += fractureVelocityMeanInCutInterpolated [ idofs[k] ];
                        }
                        else
                        {
                            fractureVelocityMeanUNCUTInterpolated [ idofs[k] ] += fractureVelocityMeanOutCutInterpolated [ idofs[k] ];
                        }

                    }

                }

            }
        }

        // fracturePressureMeanIntersectionInterpolated

		std::cout<<std::endl;
		std::cout<<std::endl;
		std::cout << "Frattura " << f <<" Pressione " << fracturePressureMeanUNCUTInterpolated << std::endl;
		std::cout<<std::endl;


		osFileName.str("");
        osFileName << "fracturePressure" << f << ".vtk";

        exportSolution ( M_exporter->getFolder() + osFileName.str(), "Pressure",
                     mfproj, fracturePressureMeanUNCUTInterpolated );

        osFileName.str("");
        osFileName << "fractureVelocity" << f << ".vtk";

        exportSolution ( M_exporter->getFolder() + osFileName.str(), "Velocity",
                     mfproj_v, fractureVelocityMeanUNCUTInterpolated );

        fracture->getData().updateU( fractureVelocityMeanUNCUTInterpolated );

    }


    return;

}// solve
