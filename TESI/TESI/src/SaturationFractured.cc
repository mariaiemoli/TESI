
#include "../include/SaturationFractured.h"

/**************************************************************************/
/*  SaturationFractured.cc                                                */
/*  Classe in cui viene risolto il problema di saturazione				  */
/**************************************************************************/


SaturationFractured::SaturationFractured( const GetPot& dataFile,
		 	 	 	 	 	 	 	 	  FracturesSetPtr_Type& fractures,
		 	 	 	 	 	 	 	 	  const BCHandlerPtr_Type& bcHandler,
		 	 	 	 	 	 	 	 	  const ExporterPtr_Type& exporter,
		 	 	 	 	 	 	 	 	  const std::string time ):
		 	 	 	 	 	 	 	 	  M_section ( time ),
		 	 	 	 	 	 	 	 	  M_fractures( fractures ),
		 	 	 	 	 	 	 	 	  M_bcHandler( bcHandler ),
		 	 	 	 	 	 	 	 	  M_exporter ( exporter ),
		 	 	 	 	 	 	 	 	  M_t ( dataFile ( ( M_section + "endTime" ).data (), 1. ) ),
		 	 	 	 	 	 	 	 	  M_fractureSaturation ( M_fractures->getNumberFractures() )
{
}

void SaturationFractured::init()
{
	size_type numberFracture = M_fractures->getNumberFractures();

    sizeVector_Type fractureNumberDOF ( numberFracture );

    // numero complessivo dei gradi di libertà
    size_type fractureTotalNumberDOF = 0;

    for ( size_type f = 0; f < numberFracture; f++ )
    {
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization() - 1;

    	fractureTotalNumberDOF += fractureNumberDOF [ f ];
    }

    // Inizializziamo tutte le matrici a blocchi, la matrice globale e il termine di destra per il sistema e il vettore delle soluzioni

    // Allochiamo la matrice globale del sistema:  M_darcyGlobalMatrix
    M_globalMatrix.reset( new sparseMatrix_Type( fractureTotalNumberDOF, fractureTotalNumberDOF ));
    gmm::clear( *M_globalMatrix );

    // Allochiamo il vettore del termine noto di destra del sistema: M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset( new scalarVector_Type( fractureTotalNumberDOF ));
    gmm::clear( *M_globalRightHandSide );

    // Allochiamo il vettore globale del termine incognito del sistema: M_darcyVelocityAndPressure
    M_saturation.reset( new scalarVector_Type( fractureTotalNumberDOF ));
    gmm::clear( *M_saturation );

    // Inizio a riempire la matrice globale del sistema
    for ( size_type i = 0; i< fractureTotalNumberDOF; i++ )
    {
    	( *M_globalMatrix )( i, i ) = 1.;
    }

    return;
}// init



void SaturationFractured::assembly( const scalar_type& landa, const scalarVectorContainer_Type& u0, const scalarVectorContainer_Type& flux )
{
	size_type numberFracture = M_fractures->getNumberFractures();

	size_type fractureShift = 0;

    sizeVector_Type fractureNumberDOF ( numberFracture );

    // numero complessivo dei gradi di libertà
    size_type fractureTotalNumberDOF = 0;

    for ( size_type f = 0; f < numberFracture; f++ )
    {
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization() -1;

    	fractureTotalNumberDOF += fractureNumberDOF [ f ];
    }

	for ( size_type f = 0; f < numberFracture; f++ )
	{
		// Impongo i flussi
		for ( size_type i = 1; i < fractureNumberDOF [ f ]; i++ )
		{
			( *M_globalRightHandSide ) [ i + fractureShift ]  = u0 [ f ][ i ] - landa*(flux [ f ] [ i ] - flux [ f ] [ i - 1 ]);
		}

		// Impongo le condizioni al contorno
        if ( M_fractures->getFracture ( f )->getData().getBc() =="f")
        {
        	scalar_type a = M_fractures->getFracture ( f )->getData().getA();
        	scalar_type b = M_fractures->getFracture ( f )->getData().getB();

        	( *M_globalRightHandSide ) [ fractureShift ]  = M_fractures->getFracture ( f )->getData().getCi( a );

        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ] = M_fractures->getFracture ( f )->getData().getCi( b );

        }
        else if ( M_fractures->getFracture ( f )->getData().getBc() =="p")
        {
        	( *M_globalRightHandSide ) [ fractureShift ]  = ( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 2 ];
        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ]  = ( *M_globalRightHandSide ) [ fractureShift + 1 ];
        }
        else
        {
        	std::cout << " bc non è stato impostato " << std::endl;
        }

        fractureShift += fractureNumberDOF [ f ];

	}

	return;

}// assembly



void SaturationFractured::solve()
{
	size_type numberFracture = M_fractures->getNumberFractures();

	scalar_type h = M_fractures->getFracture ( 0 )-> getData().getH();

	scalarVectorContainer_Type u0 ( numberFracture );

	sizeVector_Type fractureNumberDOF ( numberFracture );

	size_type fractureShift = 0;

    // numero complessivo dei gradi di libertà
    size_type fractureTotalNumberDOF = 0;

    for ( size_type f = 0; f < numberFracture; f++ )
    {
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization() - 1;

    	fractureTotalNumberDOF += fractureNumberDOF [ f ];

    }

	// Ricavo dt come il massimo valore che verifichi la condizione cfl per ogni frattura
	for ( size_type f = 0; f < numberFracture; f++ )
	{
		if ( h > M_fractures->getFracture ( f )-> getData().getH() )
		{
			h = M_fractures->getFracture ( f )-> getData().getH();
		}

		// Condizione iniziale
		u0 [ f ] = M_fractures->getFracture ( f )->getCi();

		M_fractureSaturation [ f ].reset(new scalarVector_Type( fractureNumberDOF [ f ], 0));

	}

	scalar_type dt = h/2;
	scalar_type nt = ceil(M_t/dt);
	dt = M_t/nt;
	scalar_type landa = dt/h;

	scalarVectorContainer_Type flux( numberFracture );


	std::cout << std::endl << "Solving problem in Omega..." << std::flush;

	for( size_type k = 0; k < nt; k++ )
	{
		fractureShift = 0;

		for ( size_type f = 0; f < numberFracture; f++ )
		{

			if( M_fractures->getFracture ( f )->getData().getNumberFlux() == 1 )
			{
				solve_continuity ( f, u0 [ f ], flux [ f ] );
			}
			else
			{
				solve_discontinuity ( f, u0 [ f ], flux [ f ] );
			}
		}

		assembly ( landa, u0, flux );

		// Solve the Darcy problem

		scalar_type roundConditionNumber;

		SuperLU_solve(*M_globalMatrix, *M_saturation,
					  *M_globalRightHandSide, roundConditionNumber);

		// aggiorno u0
		for (size_type f = 0; f < numberFracture; f++ )
		{
			gmm::copy( gmm::sub_vector( *M_saturation, gmm::sub_interval( fractureShift, fractureShift + fractureNumberDOF [ f ] )), u0 [ f ] );

			fractureShift += fractureNumberDOF [ f ];
		}

	}

    M_exporter->spy(M_globalMatrix, "./matlab/matrice.mm");

	fractureShift = 0;

	// esporto la soluzione
	for ( size_type f = 0; f < numberFracture; f++ )
	{
		gmm::copy ( gmm::sub_vector ( *M_saturation, gmm::sub_interval( fractureShift, fractureNumberDOF [ f ] )), *( M_fractureSaturation [ f ] ) );

		fractureShift += fractureNumberDOF [ f ];


		std::ostringstream ss;
		ss << "./matlab/saturation_" << f << ".txt";

		std::string name = ss.str();

		std::ofstream exp( name );

		exp << *( M_fractureSaturation [ f ] );


		//M_exporter->spy( M_fractureSaturation [ f ] , ss.str());

	}

	std::cout << std::endl << " completed! " << std::endl;

	return;

}// solve


void SaturationFractured::solve_continuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux )
{

	scalarVector_Type fl;
	scalarVector_Type fl1;


	// calcolo il flusso
	M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
	M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

	size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization() - 1;


	Flux.clear();
	Flux.resize( n );

	for ( size_type i = 0; i < n-1; i++)
	{
		scalar_type s = (fl[ i+1 ] - fl[ i ] )*(u0 [ i+1 ] - u0 [ i ] );

		if ( s >= 0)
		{
			Flux[ i ] = fl[ i ];
		}
		else
		{
			Flux[ i ] = fl[ i+1 ];
		}

		// entropy fix: corregge il flusso se c'è una rarefazione transonica
		if ( fl1[ i ] < 0 && fl1[ i+1] > 0)
		{
			scalar_type us;

			// trova il valore di u per il quale f'(u)=0
			if ( u0 [ i ] >= u0 [ i+1 ] )
			{
				us = M_fractures->getFracture ( f )->getData().fzero
													( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0 [ i+1 ], u0 [ i ] );
			}
			else
			{
				us = M_fractures->getFracture ( f )->getData().fzero
													( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0 [ i ], u0 [ i+1 ] );

			}

			Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( us );
		}

	}

	return;
}// solve_continuity


void SaturationFractured::solve_discontinuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux )
{
	scalar_type x_d = M_fractures->getFracture ( f )->getData().getXd();

	scalarVector_Type fl;
	scalarVector_Type fl1;

	scalarVector_Type gl;
	scalarVector_Type gl1;


	// calcolo i flussi
	M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
	M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

	M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 1 );
	M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 1 );

	size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization()-1;


	Flux.clear();
	Flux.resize( n );

	for ( size_type i = 0; i < n-1; i++)
	{
		bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_fractures->getFracture ( f )->getMeshFlat().points_of_convex ( i );

		scalar_type x1 = std::min(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);
		scalar_type x2 = std::max(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);

		if( x1 < x_d && gmm::abs(x_d - x2) > 1.0E-2 )
		{
			// siamo prima della discontinuità, usiamo il flusso 1

			scalar_type s = (fl[ i+1 ] - fl[ i ] )*(u0 [ i+1 ] - u0 [ i ] );

			if ( s >= 0)
			{
				Flux[ i ] = fl[ i ];
			}
			else
			{
				Flux[ i ] = fl[ i+1 ];
			}

			// entropy fix: corregge il flusso se c'è una rarefazione transonica
			if ( fl1[ i ] < 0 && fl1[ i+1] > 0)
			{
				scalar_type us;

				// trova il valore di u per il quale f'(u)=0
				if ( u0 [ i ] >= u0 [ i+1 ] )
				{
					us = M_fractures->getFracture ( f )->getData().fzero
														( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0 [ i+1 ], u0 [ i ] );
				}
				else
				{
					us = M_fractures->getFracture ( f )->getData().fzero
														( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0 [ i ], u0 [ i+1 ] );

				}

				Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( us );
			}

		}
		else if( x1 < x_d && gmm::abs(x_d - x2) < 1.0E-3 )
		{
			// calcolo il flusso all'interfaccia

			scalar_type s = (gl[ i+1 ] - fl[ i ] )*(u0 [ i+1 ] - u0 [ i ] );

			if ( s >= 0)
			{
				Flux[ i ] = fl[ i ];
			}
			else
			{
				Flux[ i ] = gl[ i + 1 ];
			}

		}
		else //if ( x1 > x_d /*&& gmm::abs(x_d - x1) > 1.0E-3*/ )
		{
			// siamo dopo la discontinuità, usiamo il flusso 2


			scalar_type s = (gl[ i+1 ] - gl[ i ] )*(u0 [ i+1 ] - u0 [ i ] );

			if ( s >= 0)
			{
				Flux[ i ] = gl[ i ];
			}
			else
			{
				Flux[ i ] = gl[ i+1 ];
			}

			// entropy fix: corregge il flusso se c'è una rarefazione transonica
			if ( gl1[ i ] < 0 && gl1[ i+1] > 0)
			{
				scalar_type us;

				// trova il valore di u per il quale f'(u)=0
				if ( u0 [ i ] >= u0 [ i+1 ] )
				{
					us = M_fractures->getFracture ( f )->getData().fzero
														( M_fractures->getFracture ( f )->getData().getFlux1( 1 ), u0 [ i+1 ], u0 [ i ] );
				}
				else
				{
					us = M_fractures->getFracture ( f )->getData().fzero
														( M_fractures->getFracture ( f )->getData().getFlux1( 1 ), u0 [ i ], u0 [ i+1 ] );
				}

				Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( us );
			}

		}


	}

	return;
}// solve_discontinuity
