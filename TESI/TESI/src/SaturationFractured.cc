
#include "../include/SaturationFractured.h"

/**************************************************************************/
/*  SaturationFractured.cc                                                */
/*  Classe in cui viene risolto il problema di saturazione				  */
/**************************************************************************/


SaturationFractured::SaturationFractured( const GetPot& dataFile,
		 	 	 	 	 	 	 	 	  FracturesSetPtr_Type& fractures,
		 	 	 	 	 	 	 	 	  const std::string time ):
		 	 	 	 	 	 	 	 	  M_section ( time ),
		 	 	 	 	 	 	 	 	  M_t ( dataFile ( ( M_section + "endTime" ).data (), 1. ) ),
		 	 	 	 	 	 	 	 	  M_fractures( fractures ),
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
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization();

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
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization();

    	fractureTotalNumberDOF += fractureNumberDOF [ f ];
    }


	for ( size_type f = 0; f < numberFracture; f++ )
	{
		// Impongo i flussi
		for ( size_type i = 1; i < fractureNumberDOF [ f ] - 1; i++ )
		{
			( *M_globalRightHandSide ) [ i + fractureShift ]  = u0 [ f ][ i ] - landa*(flux [ f ] [ i ] - flux [ f ] [ i - 1 ]);
		}

		// Impongo le condizioni al contorno
        if ( M_fractures->getFracture ( f )->getData().getBc() =="f")
        {
        	( *M_globalRightHandSide ) [ fractureShift ]  = M_fractures->getFracture ( f )->getData().getUl();

        	( *M_globalRightHandSide ) [ fractureShift + fractureNumberDOF [ f ] - 1 ] = M_fractures->getFracture ( f )->getData().getUr();

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

	std::cout << " numberFracture:  " << numberFracture << std::endl;

	sizeVector_Type fractureNumberDOF ( numberFracture );

	size_type fractureShift = 0;

    // numero complessivo dei gradi di libertà
    size_type fractureTotalNumberDOF = 0;

    for ( size_type f = 0; f < numberFracture; f++ )
    {
    	fractureNumberDOF [ f ] = M_fractures ->getFracture ( f )->getData().getSpatialDiscretization();

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

	scalarVector_Type fl;
	scalarVector_Type fl1;
	scalarVectorContainer_Type flux( numberFracture );


	for( size_type k = 0; k < nt; k++ )
	{

		for ( size_type f = 0; f < numberFracture; f++ )
		{

			// calcolo il flusso
			M_fractures->getFracture ( f )->getData().feval( fl, u0 [ f ], 1 );

			size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization();


			flux [ f ].clear();
			flux [ f ].resize( n );

			for ( size_type i = 0; i < n-2; i++)
			{
		        scalar_type s = (fl[ i+1 ] - fl[ i ] )*(u0 [ f ][ i+1 ] - u0 [ f ][ i ] );

		        if ( s >= 0)
		        {
		        	flux [ f ][ i ] = fl[ i ];
		        }
		        else
		        {
		        	flux [ f ][ i ] = fl[ i+1 ];
		        }

		        // entropy fix: corregge il flusso se c'e' una rarefazione transonica

		        M_fractures->getFracture ( f )->getData().feval( fl1, u0 [ f ], 2 );

		        if ( fl1[ i ] < 0 && fl1[ i+1] > 0)
		        {
		        	scalar_type us;

		        	// trova il valore di u per il quale f'(u)=0
		        	if ( u0 [ f ][ i ] >= u0 [ f ][ i+1 ] )
		        	{
		        		us = M_fractures->getFracture ( f )->getData().fzero
		        											( M_fractures->getFracture ( f )->getData().getFlux1(), u0 [ f ][ i+1 ], u0 [ f ][ i ] );
		        	}
		        	else
		        	{
		        		us = M_fractures->getFracture ( f )->getData().fzero
		        											( M_fractures->getFracture ( f )->getData().getFlux1(), u0 [ f ][ i ], u0 [ f ][ i+1 ] );

		        	}

		        	flux [ f ][ i ] = M_fractures->getFracture ( f )->getData().feval_scal( us );
		        }

			}


		}

		assembly ( landa, u0, flux );

		// risolvo
	    // Solve the Darcy problem
	    std::cout << std::endl << "Solving problem in Omega..." << std::flush;
	    scalar_type roundConditionNumber;

	    SuperLU_solve(*M_globalMatrix, *M_saturation,
	                  *M_globalRightHandSide, roundConditionNumber);


		// aggiorno u0
        //gmm::copy( gmm::sub_vector( *M_saturation, gmm::sub_interval( 0, fractureTotalNumberDOF )), u0 );


	    for (size_type f = 0; f < numberFracture; f++ )
	    {
	    	gmm::copy( gmm::sub_vector( *M_saturation, gmm::sub_interval( fractureShift, fractureShift + fractureNumberDOF [ f ] )), u0 [ f ] );
	    }

	}

	fractureShift = 0;

	// esporto la soluzione
	for ( size_type f = 0; f < numberFracture; f++ )
	{
		gmm::copy ( gmm::sub_vector ( *M_saturation, gmm::sub_interval( fractureShift, fractureNumberDOF [ f ] )), *( M_fractureSaturation [ f ] ) );

		fractureShift += fractureNumberDOF [ f ];

		std::ostringstream ss;
		ss << "saturation_" << f << ".txt";
		std::string name = ss.str();

		std::ofstream exp( name );

		exp << *( M_fractureSaturation [ f ] );

	}


	//	std::ostringstream ss;
	//	ss << i;
	//	std::string name("u" );
	//	name = name + ss.str();
	//	std::ofstream exp("u.txt");

	//	exp << u;




	return;
}// solve


