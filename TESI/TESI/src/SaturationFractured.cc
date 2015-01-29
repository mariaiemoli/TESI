
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

	//size_type NumIntersection = M_fractures->getIntersections()->getNumberIntersection ();
    size_type NumIntersection = 0;

    // Inizializziamo tutte le matrici a blocchi, la matrice globale e il termine di destra per il sistema e il vettore delle soluzioni

    // Allochiamo la matrice globale del sistema:  M_darcyGlobalMatrix
    M_globalMatrix.reset( new sparseMatrix_Type( fractureTotalNumberDOF + NumIntersection, fractureTotalNumberDOF + NumIntersection ));
    gmm::clear( *M_globalMatrix );

    // Allochiamo il vettore del termine noto di destra del sistema: M_darcyGlobalRightHandSide
    M_globalRightHandSide.reset( new scalarVector_Type( fractureTotalNumberDOF + NumIntersection ));
    gmm::clear( *M_globalRightHandSide );

    // Allochiamo il vettore globale del termine incognito del sistema: M_darcyVelocityAndPressure
    M_saturation.reset( new scalarVector_Type( fractureTotalNumberDOF + NumIntersection ));
    gmm::clear( *M_saturation );

    // Inizio a riempire la matrice globale del sistema
    for ( size_type i = 0; i< fractureTotalNumberDOF + NumIntersection; i++ )
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
		/*
		scalar_type flux_in;
		scalar_type flux_out;

		size_type H = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getH();
		scalar_type Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();

		std::string monotone = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getMonotone();

		if ( monotone == "false" )
		{
			if ( H == 2. )
			{
				scalar_type m1 = std::max( u0[ f ][ 0 ], Us );
				scalar_type m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUin();

				scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m2, 0 );
				scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );

				flux_in = std::max( F_ul, F_ur );

				m1 = std::max( u0[ f ][ fractureNumberDOF [ f ] - 1 ], Us );
				m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUout();


				F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, 0 );


				flux_out = std::max( F_ul, F_ur );

			}
			else
			{
				scalar_type m1 = std::min( u0[ f ][ 0 ], Us );
				scalar_type m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUin();

				scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m2, 0 );
				scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );

				flux_in = std::min( F_ul, F_ur );

				m1 = std::min( u0[ f ][ fractureNumberDOF [ f ] - 1 ], Us );
				m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUout();

				F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, 0 );

				flux_out = std::min( F_ul, F_ur );

			}
		}
		else
		{
			scalar_type Sl = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUin();

			std::cout << " Sl " << Sl << std::endl;

			scalar_type Sr;

			Sr = u0[ f ][ 0 ];

			scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, 0 );
			scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, 0 );

			flux_in = F_ul;
		}
		*/
		//( *M_globalRightHandSide ) [ fractureShift ]  = u0 [ f ][ 0 ] - landa*(flux [ f ] [ 0 ] - flux_in );

		//( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ]-1 ]  = u0 [ f ][ fractureNumberDOF [ f ] -1 ] - landa*(flux_out - flux [ f ] [ fractureNumberDOF [ f ] - 2 ] );

		//std::cout << "  frattura  " << f << "  flux_in  " << flux_in << "   flux_out " << flux_out << std::endl;


		int U = M_fractures->getFracture ( f )->getData().getVelocity();

	//	if ( U > 0 )
		{
			// constrollare qui
			// Impongo i flussi
			for ( size_type i = 0; i < fractureNumberDOF [ f ] ; i++ )
			{
				( *M_globalRightHandSide ) [ i + fractureShift ]  = u0 [ f ][ i ] - landa*(flux [ f ] [ i+1 ] - flux [ f ] [ i  ]);
			}
		}
	/*	else
		{
			// Impongo i flussi
			for ( size_type i = 0; i < fractureNumberDOF [ f ] ; i++ )
			{
				( *M_globalRightHandSide ) [ i + fractureShift ]  = u0 [ f ][ i ] - landa*(flux [ f ] [ i+1 ] - flux [ f ] [ i  ]);
			}
		}*/
		/*
    	scalar_type a = M_fractures->getFracture ( f )->getData().getA();
    	scalar_type b = M_fractures->getFracture ( f )->getData().getB();


		// Impongo le condizioni al contorno: può essere che io voglia forzare il valore di saturazione nel primo e nell'ultimo nodo
        if ( M_fractures->getFracture ( f )->getData().getBc() =="f")
        {
        	( *M_globalRightHandSide ) [ fractureShift ]  = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUin();

        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ] = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUout();


        }
        else if ( M_fractures->getFracture ( f )->getData().getBc() =="p")
        {
        	( *M_globalRightHandSide ) [ fractureShift ]  = ( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 2 ];
        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ]  = ( *M_globalRightHandSide ) [ fractureShift + 1 ];
        }
		*/
        fractureShift += fractureNumberDOF [ f ];

	}


	return;

}// assembly



void SaturationFractured::solve()
{
	size_type numberFracture = M_fractures->getNumberFractures();

	scalar_type h = M_fractures->getFracture ( 0 )-> getData().getH();
	scalar_type cfl = M_fractures->getFracture ( 0 )->getData().getCfl();

	scalarVectorContainer_Type u0 ( numberFracture );

	sizeVector_Type fractureNumberDOF ( numberFracture );

	size_type fractureShift = 0;

	size_type NumIntersection = M_fractures->getIntersections()->getNumberIntersection ();

	IntersectDataContainer_Type Intersection = M_fractures->getIntersections()-> getIntersections ();

	scalar_type dt = h/cfl;
	scalar_type nt = ceil(M_t/dt);
	dt = M_t/nt;
	scalar_type landa = dt/h;


	nt = 50;
	dt = 0.002;

	landa = dt/h;

    // risolvo per ogni intersezione
	for ( size_type i = 0; i < NumIntersection; i++ )
	{
		/*
		 * Per ogni intersezione calcolo u_I, così facendo aggiorno anche i di  Uin e valori per ogni frattura
		 */
		Intersection [ i ].update_Ui( landa );

		std::cout << " M_u0 " << Intersection [ i ].getU0() << std::endl;

	}

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

		if ( cfl < M_fractures->getFracture ( f )->getData().getCfl() )
		{
			cfl = M_fractures->getFracture ( f )->getData().getCfl();
		}

		// Condizione iniziale
		u0 [ f ] = M_fractures->getFracture ( f )->getCi();

		M_fractureSaturation [ f ].reset(new scalarVector_Type( fractureNumberDOF [ f ], 0));

	}

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

			std::ostringstream ss;
			ss << "./matlab/risultati/saturation_flusso.txt";

			std::string name = ss.str();

			std::ofstream exp( name );

			exp << flux[ f ];
		}

		assembly ( landa, u0, flux );

		// Solve the Saturation problem
		scalar_type roundConditionNumber;

		SuperLU_solve(*M_globalMatrix, *M_saturation,
					  *M_globalRightHandSide, roundConditionNumber);


		// aggiorno u0
		for (size_type f = 0; f < numberFracture; f++ )
		{
			gmm::copy( gmm::sub_vector( *M_saturation, gmm::sub_interval( fractureShift, fractureNumberDOF [ f ] )), u0 [ f ] );

			fractureShift += fractureNumberDOF [ f ];
		}

	}

    M_exporter->spy( M_globalMatrix, "./matlab/matrice.mm" );

	fractureShift = 0;


	// esporto la soluzione
	for ( size_type f = 0; f < numberFracture; f++ )
	{
		gmm::copy ( gmm::sub_vector ( *M_saturation, gmm::sub_interval( fractureShift, fractureNumberDOF [ f ] )), *( M_fractureSaturation [ f ] ) );

		fractureShift += fractureNumberDOF [ f ];


		std::ostringstream ss;
		ss << "./matlab/risultati/saturation_" << f << ".txt";

		std::string name = ss.str();

		std::ofstream exp( name );

		exp << *( M_fractureSaturation [ f ] );

	}

	std::cout << " completed!" << std::endl;
	std::cout << std::endl;


	return;

}// solve

void SaturationFractured::solve_continuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux, const size_type k )
{
	scalar_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization()-1;

	Flux.clear();
	Flux.resize( n );

	scalarVector_Type U0 = u0;

	std::string monotone = M_fractures->getFracture ( f )->getData().getFluxHandler( k )->getMonotone();

	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type Us = 0.;

	size_type H = M_fractures->getFracture ( f )->getData().getFluxHandler( k )->getH();

	if( monotone == "false" )
	{
		Us = M_fractures->getFracture ( f )->getData().getFluxHandler( k )->getUs();
	}


	for ( size_type i = 0; i < n; i++)
	{
		if ( monotone == "false" )
		{
			if( H == 2.0 )		// caso A
			{
				if ( i != n-1 )
				{
					scalar_type m1 = std::max( u0[ i ], Us );
					scalar_type m2 = std::min( u0[ i+1 ], Us );

					scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, k );
					scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

					Flux[ i ] = std::max( F_ul, F_ur );
				}
				else
				{
					scalar_type m1 = std::max( u0[ i ], Us );
					scalar_type m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( k )->getUout();

					scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, k );
					scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

					Flux[ i ] = std::max( F_ul, F_ur );
				}
			}
			else		// caso B
			{
				if ( i != n-1 )
				{
					scalar_type m1 = std::min( u0[ i ], Us );
					scalar_type m2 = std::max( u0[ i+1 ], Us );

					scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, k );
					scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

					Flux[ i ] = std::min( F_ul, F_ur );
				}
				else
				{
					scalar_type m1 = std::min( u0[ i ], Us );
					scalar_type m2 = M_fractures->getFracture ( f )->getData().getFluxHandler( k )->getUout();

					scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( m1, k );
					scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

					Flux[ i ] = std::min( F_ul, F_ur );

				}
			}
		}
		else
		{
			// implemento il metodo di anna

			scalar_type Sl = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUin();

			int U = M_fractures->getFracture ( f )->getData().getVelocity();

			scalar_type Sr;


			if ( U < 0 )
			{
				U0 [ n-1 ] = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUout();
			}

			for ( size_type i = 0; i < n; i++ )
			{
				Sr = U0[ i ];

				scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, 0 );
				scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, 0 );

				if ( U > 0 )
				{
					Flux[ i ] = F_ul;
				}
				else
				{
					Flux[ i ] = F_ur;
				}

				Sl = Sr;

			}
		}
	}

	return;
}// solve_continuity


void SaturationFractured::solve_discontinuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux )
{
	size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization() - 1;
	scalar_type x_d = M_fractures->getFracture ( f )->getData().getXd();

	Flux.clear();
	Flux.resize( n );

	scalarVector_Type fl;
	scalarVector_Type fl1;

	scalarVector_Type gl;
	scalarVector_Type gl1;

	scalarVector_Type F(n);
	scalarVector_Type G(n);

	size_type H = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getH();


	// calcolo i flussi
	M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
	M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

	M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 1 );
	M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 1 );


	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type S_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
	scalar_type D_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getUs();

	// flusso 1
	solve_continuity ( f, u0, F, 0 );

	// flusso 2
	solve_continuity ( f, u0, G, 1 );

	F[ n-1 ]= fl[ n-1 ];
	G[ n-1 ]= gl[ n-1 ];

	for ( size_type i = 0; i < n-1; i++)
	{
		bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_fractures->getFracture ( f )->getMeshFlat().points_of_convex ( i );

		scalar_type x1 = std::min(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);
		scalar_type x2 = std::max(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);
		if ( x1 < x_d &&  x2 != x_d )
		{
			Flux[ i ] = F[ i ];
		}
		else if ( x1 < x_d && gmm::abs( x2-x_d ) < 1.0E-5 ) // x2 == x_d )
		{
			if ( H == 2. )
			{
				scalar_type m1 = std::max( u0[ i ], S_Us );
				scalar_type m2 = std::min( u0[ i+1 ], D_Us );

				scalar_type F_m1 = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				scalar_type F_m2 = M_fractures->getFracture ( f )->getData().feval_scal( m2, 1 );

				Flux[ i ] = std::max( F_m1, F_m2 );
			}
			else
			{
				scalar_type m1 = std::min( u0[ i ], S_Us );
				scalar_type m2 = std::max( u0[ i+1 ], D_Us );

				scalar_type F_m1 = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				scalar_type F_m2 = M_fractures->getFracture ( f )->getData().feval_scal( m2, 1 );

				Flux[ i ] = std::min( F_m1, F_m2 );
			}

		}
		else
		{
			Flux[ i ] = G[ i ];
		}


	}

	Flux[ n-1 ] = G[ n-1 ];



	return;

}// solve_discontinuity
