
#include "../include/SaturationFractured.h"

/**************************************************************************/
/*  SaturationFractured.cc                                                */
/*  Classe in cui viene risolto il problema di saturazione				  */
/**************************************************************************/


SaturationFractured::SaturationFractured( const GetPot& dataFile,
		 	 	 	 	 	 	 	 	  FracturesSetPtr_Type& fractures,
		 	 	 	 	 	 	 	 	  const ExporterPtr_Type& exporter,
		 	 	 	 	 	 	 	 	  const std::string time ):
		 	 	 	 	 	 	 	 	  M_section ( time ),
		 	 	 	 	 	 	 	 	  M_fractures( fractures ),
		 	 	 	 	 	 	 	 	  M_exporter ( exporter ),
		 	 	 	 	 	 	 	 	  M_t ( dataFile ( ( M_section + "endTime" ).data (), 1. ) ),
		 	 	 	 	 	 	 	 	  M_dt ( dataFile ( ( M_section + "dt" ).data (), 0.001 ) ),
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



void SaturationFractured::assembly( const scalarVector_Type& landa, const scalarVectorContainer_Type& u0, const scalarVectorContainer_Type& flux )
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
		for ( size_type i = 0; i < fractureNumberDOF [ f ] ; i++ )
		{
			( *M_globalRightHandSide ) [ i + fractureShift ]  = u0 [ f ][ i ] - landa[ f ]*(flux [ f ] [ i+1 ] - flux [ f ] [ i ]);
		}

        fractureShift += fractureNumberDOF [ f ];
	}

	return;

}// assembly



void SaturationFractured::solve()
{
	size_type numberFracture = M_fractures->getNumberFractures();

	scalarVectorContainer_Type u0 ( numberFracture );

	sizeVector_Type fractureNumberDOF ( numberFracture );

	size_type fractureShift = 0;

	size_type NumIntersection = M_fractures->getIntersections()->getNumberIntersection ();

	IntersectDataContainer_Type Intersection = M_fractures->getIntersections()-> getIntersections ();

	scalar_type nt = M_t/M_dt;

	scalarVector_Type landa ( numberFracture );

	for ( size_type i = 0; i < numberFracture; i++ )
	{
		scalar_type h = M_fractures->getFracture ( i )-> getH();

		landa [ i ] = M_dt/h; //0.2; //M_dt/h;
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
			solve_riemann ( f, u0 [ f ], flux [ f ] );

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

	    // risolvo per ogni intersezione
		for ( size_type i = 0; i < NumIntersection; i++ )
		{
			/*
			 * Per ogni intersezione calcolo u_I, così facendo aggiorno anche i di  Uin e valori per ogni frattura
			 */
			Intersection [ i ].update_Ui( M_dt );
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


void SaturationFractured::solve_riemann ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux )
{
	size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization() - 1;
	scalar_type x_d = M_fractures->getFracture ( f )->getData().getXd();

	Flux.clear();
	Flux.resize( n + 1 );

	scalarVector_Type fl;
	scalarVector_Type fl1;

	scalarVector_Type gl;
	scalarVector_Type gl1;

	scalarVector_Type F(n);
	scalarVector_Type G(n);

	size_type H = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getH();

	size_type numberFlux = M_fractures->getFracture ( f )->getData().getNumberFlux ();

	size_type k = 1;

	scalar_type S_Us;
	scalar_type D_Us;

	if ( numberFlux == 2 )
	{
		// calcolo i flussi
		M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
		M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

		M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 1 );
		M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 1 );


		// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
		S_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
		D_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getUs();
	}
	else
	{
		// calcolo i flussi
		M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
		M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

		M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 0 );
		M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 0 );


		// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
		S_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
		D_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();

		k = 0;
	}

	// flusso 1
	solve_continuity ( f, u0, F, 0 );

	// flusso 2
	solve_continuity ( f, u0, G, k );

	for ( size_type i = 0; i < n; i++)
	{
		bgeot::basic_mesh::ref_mesh_pt_ct nodes = M_fractures->getFracture ( f )->getMeshFlat().points_of_convex ( i );

		scalar_type x1 = std::min(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);
		scalar_type x2 = std::max(nodes [ 0 ] [ 0 ], nodes [ 1 ] [ 0 ]);

		if ( x1 < x_d &&  x2 != x_d )
		{
			Flux[ i ] = F[ i ];
		}
		else if ( x1 < x_d && gmm::abs( x2-x_d ) < 1.0E-5 )
		{
			if ( H == 2. )
			{
				scalar_type m1 = std::max( u0[ i ], S_Us );
				scalar_type m2 = std::min( u0[ i+1 ], D_Us );

				scalar_type F_m1 = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				scalar_type F_m2 = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

				Flux[ i ] = std::max( F_m1, F_m2 );
			}
			else
			{
				scalar_type m1 = std::min( u0[ i ], S_Us );
				scalar_type m2 = std::max( u0[ i+1 ], D_Us );

				scalar_type F_m1 = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
				scalar_type F_m2 = M_fractures->getFracture ( f )->getData().feval_scal( m2, k );

				Flux[ i ] = std::min( F_m1, F_m2 );
			}

		}
		else
		{
			Flux[ i ] = G[ i ];
		}

	}

	Flux[ n ] = G[ n ];

	return;

}// solve_riemann


void SaturationFractured::solve_continuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux, const size_type k )
{
	scalar_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization()-1;

	Flux.clear();
	Flux.resize( n+1 );

	std::string monotone = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getMonotone();

	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type Us = 0.;

	size_type H = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getH();

	if( monotone == "false" )
	{
		Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
	}

	if ( monotone == "false" )
	{
		scalar_type Sl;
		scalar_type Sr;

		if ( H == 2 )
		{
			Sl = std::max( M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->geta(), Us );
		}
		else
		{
			Sl = std::min( M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->geta(), Us );
		}

		for ( size_type i = 0; i < n; i++)
		{

			if( H == 2.0 )		// caso A
			{

				Sr = std::min( u0[ i ], Us );

				scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, k );
				scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, k );

				Flux[ i ] = std::max( F_ul, F_ur );

				Sl = std::max( u0[ i ], Us );
			}
			else		// caso B
			{
				Sr = std::max( u0[ i ], Us );

				scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, k );
				scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, k );

				Flux[ i ] = std::min( F_ul, F_ur );

				Sl = std::min( u0[ i ], Us );
			}
		}

		if ( H == 2 )
		{
			Sr = std::min( M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getb(), Us );
		}
		else
		{
			Sr = std::max( M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getb(), Us );
		}

		scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, k );
		scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, k );


		if ( H == 2 )
		{
			Flux[ n ] = std::max( F_ul, F_ur );
		}
		else
		{
			Flux[ n ] = std::min( F_ul, F_ur );
		}


		M_fractures->getFracture ( f )->getData().updateSI ( Flux );

	}
	else
	{
		// implemento il metodo di anna

		scalar_type Sl = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->geta();

		scalar_type Sr;

		for ( size_type i = 0; i < n; i++ )
		{
			Sr = u0[ i ];

			scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, 0 );
			scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, 0 );

			scalar_type s = ( F_ur - F_ul )*( Sr - Sl );

			if ( s > 0 )
			{
				Flux[ i ] = F_ul;
			}
			else
			{
				Flux[ i ] = F_ur;
			}

			Sl = Sr;

		}

		Sr = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getb();

		scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( Sr, k );
		scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( Sl, k );

		scalar_type s = ( F_ur - F_ul )*( Sr - Sl );

		if ( s > 0 )
		{
			Flux[ n ] = F_ul;
		}
		else
		{
			Flux[ n ] = F_ur;
		}

		M_fractures->getFracture ( f )->getData().updateSI ( Flux );
	}

	return;
}// solve_continuity
