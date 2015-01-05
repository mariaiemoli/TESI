
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

    	scalar_type a = M_fractures->getFracture ( f )->getData().getA();
    	scalar_type b = M_fractures->getFracture ( f )->getData().getB();

		// Impongo le condizioni al contorno
        if ( M_fractures->getFracture ( f )->getData().getBc() =="f")
        {
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
        	//std::cout << " bc non è stato impostato " << std::endl;

        	/*
        	// questo non va bene, va sistemato
        	if ( M_fractures->getFracture ( f )->getData().feval_scal( a, 1 ) > 0 )
        	{
            	( *M_globalRightHandSide ) [ fractureShift ]  = M_fractures->getFracture ( f )->getData().getCi( a );

        	}
        	if ( M_fractures->getFracture ( f )->getData().feval_scal( b, 1 ) > 0 )
        	{
            	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ]  = M_fractures->getFracture ( f )->getData().getCi( b );
        	}
        	*/

        	( *M_globalRightHandSide ) [ fractureShift ]  = 0.25;
        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ]  = 0.66;
        }

        fractureShift += fractureNumberDOF [ f ];

	}

	size_type NumCross = M_fractures->getIntersections()->getNumberCross ();
	size_type NumBifurcation = M_fractures->getIntersections()->getNumberBifurcation ();


	for ( size_type i = 0; i < NumBifurcation; i++ )
	{
		/*
		 * devo eliminare la riga corrispondente della riga e l'elemento corrispondente del termine noto per ogni frattura
		 */

		std::cout << " ciao " << std::endl;
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



	scalar_type dt = h/cfl;
	scalar_type nt = ceil(M_t/dt);
	dt = M_t/nt;
	scalar_type landa = dt/h;


	/*
	scalar_type dt = 0.02;
	scalar_type nt = ceil(M_t/dt);
	scalar_type landa = dt/h;
	*/

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

    M_exporter->spy(M_globalMatrix, "./matlab/matrice.mm");

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


void SaturationFractured::solve_continuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux, const size_type i )
{
	size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization() - 1;

	Flux.clear();
	Flux.resize( n );

	std::string monotone = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getMonotone();

	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type Us;

	scalar_type F_Us;

	if( monotone == "false" )
	{
		Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();

		F_Us = M_fractures->getFracture ( f )->getData().feval_scal( Us, 0 );
	}


	for ( size_type i = 0; i < n-1; i++)
	{
		if ( monotone == "false" )
		{
			if( F_Us >= 0 )		// caso A
			{
				if ( u0[ i ] <= Us )		// caso 1
				{
					if ( u0[ i+1 ] <= Us )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
					}
					else
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 0 );
						Flux[ i ] = std::min( F_ul, F_ur );
					}
				}
				else					// caso 2
				{
					if ( u0[ i+1 ] <= Us )
					{
						Flux[ i ] = F_Us;
					}
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 0 );
					}
				}
			}
			else		// caso B
			{
				if ( u0[ i ] <= Us )		// caso 1
				{
					if ( u0[ i+1 ] <= Us )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 0 );
					}
					else
					{
						Flux[ i ] = F_Us;
					}
				}
				else				// caso 2
				{
					if ( u0[ i+1 ] <= Us )
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 0 );
						Flux[ i ] = std::max( F_ul, F_ur );
					}
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
					}
				}

			}
		}
		else
		{
			scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
			scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 0 );

			scalar_type s = (F_ur - F_ul )*(u0[ i+1 ] - u0[ i ] );


	        if ( s >= 0)
	        {
	        	Flux[ i ] = F_ul;
	        }
	        else
	        {
	        	Flux[ i ] = F_ur;
	        }

	        // entropy fix: corregge il flusso se c'e' una rarefazione transonica

	        scalarVector_Type fl1;

	        M_fractures->getFracture ( f )->getData().feval( fl1, u0, 0, 0 );

	        if ( fl1[ i ] < 0 && fl1[ i+1 ] > 0)
	        {
	        	scalar_type us;

	        	// trova il valore di u per il quale f'(u)=0
	        	if ( u0[ i ] >= u0[ i+1 ] )
	        	{
	        		us = M_fractures->getFracture ( f )->getData().fzero ( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0[ i+1 ], u0[ i ] );
	        	}
	        	else
	        	{
	        		us = M_fractures->getFracture ( f )->getData().fzero ( M_fractures->getFracture ( f )->getData().getFlux1( 0 ), u0[ i ], u0[ i+1 ] );

	        	}

	        	Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( us, 0, 0 );
	        }
		}
	}

	return;
}// solve_continuity

/*
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


	// calcolo i flussi
	M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
	M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

	M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 1 );
	M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 1 );


	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type S_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
	scalar_type D_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getUs();

	scalar_type SF_Us = M_fractures->getFracture ( f )->getData().feval_scal( S_Us, 0 );
	scalar_type DF_Us = M_fractures->getFracture ( f )->getData().feval_scal( D_Us, 1 );

	// flusso 1
	solve_continuity ( f, u0, F );

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
			scalar_type m1 = std::max( u0[ i ], S_Us );
			scalar_type m2 = std::min( u0[ i+1 ], D_Us );

			scalar_type F_m1 = M_fractures->getFracture ( f )->getData().feval_scal( m1, 0 );
			scalar_type F_m2 = M_fractures->getFracture ( f )->getData().feval_scal( m2, 1 );

			Flux[ i ] = std::max( F_m1, F_m2 );

		}
		else
		{
			Flux[ i ] = G[ i ];
		}


	}

	Flux[ n-1 ] = G[ n-1 ];



	return;

}// solve_discontinuity
*/


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


	// calcolo i flussi
	M_fractures->getFracture ( f )->getData().feval( fl, u0, 1, 0 );
	M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2, 0 );

	M_fractures->getFracture ( f )->getData().feval( gl, u0, 1, 1 );
	M_fractures->getFracture ( f )->getData().feval( gl1, u0, 2, 1 );


	// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
	scalar_type S_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getUs();
	scalar_type D_Us = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getUs();

	scalar_type S_first = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getFirst();
	scalar_type D_first = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getFirst();

	scalar_type S_second = M_fractures->getFracture ( f )->getData().getFluxHandler( 0 )->getSecond();
	scalar_type D_second = M_fractures->getFracture ( f )->getData().getFluxHandler( 1 )->getSecond();

	scalar_type SF_Us = M_fractures->getFracture ( f )->getData().feval_scal( S_Us, 0 );
	scalar_type DF_Us = M_fractures->getFracture ( f )->getData().feval_scal( D_Us, 1 );

	// flusso 1
	solve_continuity ( f, u0, F );

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
		else if ( x1 < x_d && gmm::abs( x2-x_d ) < 1.0E-5 )  // x2 == x_d )
		{
			// caso 1
			if ( SF_Us <= 0 && SF_Us >= DF_Us )
			{
				// caso 1.1
				if ( u0[ i ] < S_Us )
				{
					// caso 1.1.1
					if ( u0[ i+1 ] < D_first )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}
					// caso 1.1.2
					else
					{
						Flux[ i ] = SF_Us;
					}
				}
				// caso 1.2
				else
				{
					// caso 1.2.1
					if ( u0[ i+1 ] < D_first )
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
						Flux[ i ] = std::max( F_ul, F_ur );
					}
					// caso 1.2.2
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
					}

				}
			} // caso 1

			// caso 2
			else if ( DF_Us < 0 && SF_Us < DF_Us )
			{
				// caso 2.1
				if ( u0[ i ] < S_second )
				{
					// caso 2.1.1
					if ( u0[ i+1 ] < D_Us )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}
					// caso 2.1.2
					else
					{
						Flux[ i ] = DF_Us;
					}
				}
				// caso 2.2
				else
				{
					// caso 2.2.1
					if ( u0[ i+1 ] < D_Us )
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
						Flux[ i ] = std::max( F_ul, F_ur );
					}
					// caso 2.2.2
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
					}

				}
			} // caso 2

			// caso 3
			else if ( SF_Us > 0 && SF_Us < DF_Us )
			{
				// caso 3.1
				if ( u0[ i ] < S_Us )
				{
					// caso 3.1.1
					if ( u0[ i+1 ] < D_second )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}
					// caso 3.1.2
					else
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 0 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
						Flux[ i ] = std::min( F_ul, F_ur );
					}
				}
				// caso 3.2
				else
				{
					// caso 3.2.1
					if ( u0[ i+1 ] < D_second )
					{
						Flux[ i ] = DF_Us;
					}
					// caso 3.2.2
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}

				}
			} // caso 3

			// caso 4
			else if ( DF_Us > 0 && SF_Us > DF_Us )
			{
				// caso 4.1
				if ( u0[ i ] < S_first )
				{
					// caso 4.1.1
					if ( u0[ i+1 ] < D_Us )
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}
					// caso 4.1.2
					else
					{
						scalar_type F_ul = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i ], 1 );
						scalar_type F_ur = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
						Flux[ i ] = std::min( F_ul, F_ur );
					}
				}
				// caso 4.2
				else
				{
					// caso 4.2.1
					if ( u0[ i+1 ] < D_Us )
					{
						Flux[ i ] = DF_Us;
					}
					// caso 4.2.2
					else
					{
						Flux[ i ] = M_fractures->getFracture ( f )->getData().feval_scal( u0[ i+1 ], 1 );
					}

				}
			} // caso 4

		}
		else
		{
			Flux[ i ] = - G[ i ];
		}


	}

	Flux[ n-1 ] = G[ n-1 ];



	return;

}// solve_discontinuity
