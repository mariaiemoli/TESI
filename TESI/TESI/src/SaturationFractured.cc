
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
		 	 	 	 	 	 	 	 	  M_fractures( fractures )
{
}

void SaturationFractured::solve()
{
	size_type numberFracture = M_fractures->getNumberFractures();

	for ( size_type f = 0; f < numberFracture; f++ )
	{
		size_type n = M_fractures->getFracture ( f )-> getData().getSpatialDiscretization();

		// passo della mesh
		scalar_type h = M_fractures->getFracture ( f )-> getData().getH();

		scalar_type landa0=0.9;
		scalar_type cfl = M_fractures->getFracture ( f )->getData().getCfl();
		scalar_type dt=landa0*h/(cfl);

		scalar_type nt = ceil(M_t/dt);
		dt = M_t/nt;
		scalar_type landa = dt/h;

		std::cout << " dt:  " << dt << " cfl:  " << cfl << " h:  " << h << std::endl;
		std::cout << " dt*cfl:  " << dt*cfl << "   h/2:  " << 0.5*h << std::endl;
		std::cout << "landa:  " << landa << std::endl;

		scalarVector_Type fl;
		scalarVector_Type fl1;
		scalarVector_Type flux;

		// ciclo temporale
		scalarVector_Type u0 = M_fractures->getFracture ( f )->getCi();

		std::cout << "u0.size() " << u0.size() << std::endl;
		scalarVector_Type u, uu0;

		uu0 = u0;

		std::cout << "  NT:  " << nt << std::endl;

		for ( size_type i = 0; i < nt; i++)
		{
			std::cout << " passo " << i << std::endl;
			// calcolo il flusso: memorizza f(i+1/2) in f(i)
			M_fractures->getFracture ( f )->getData().feval( fl, u0, 1 );
			flux.clear();
			flux.resize( n );

			u.clear();
			u.resize( n );

			for ( size_type j = 0; j < n-2; j++)
			{
		        scalar_type s = (fl[ j+1 ] - fl[ j] )*(u0[ j+1 ] - u0[ j] );

		        if ( s >= 0)
		        {
		        	flux[ j ] = fl[ j ];
		        }
		        else
		        {
		        	flux[ j ] = fl[ j+1 ];
		        }

		        // entropy fix: corregge il flusso se c'e' una rarefazione transonica

		        M_fractures->getFracture ( f )->getData().feval( fl1, u0, 2 );

		        if ( fl1[ j ] < 0 && fl1[ j+1] > 0)
		        {
		        	scalar_type us;

		        	// trova il valore di u per il quale f'(u)=0
		        	if ( u0[ j ] >= u0[ j+1 ] )
		        	{
		        		us = M_fractures->getFracture ( f )->getData().fzero ( M_fractures->getFracture ( f )->getData().getFlux1(), u0[ j+1 ], u0[ j ] );
		        	}
		        	else
		        	{
		        		us = M_fractures->getFracture ( f )->getData().fzero ( M_fractures->getFracture ( f )->getData().getFlux1(), u0[ j ], u0[ j+1 ] );

		        	}

		        	flux[ j ] = M_fractures->getFracture ( f )->getData().feval_scal( us );
		        }

			}


	        // calcolo la soluzione
	        for ( size_type j = 1; j < n-1; j++ )
	        {
	        	u[ j ] = u0[ j ] - landa*(flux[ j ]-flux[ j-1 ]);

	        }

	        // calcolo le condizioni al contorno
	        if ( M_fractures->getFracture ( f )->getData().getBc() =="f")
	        {
	        		u[ 0 ] = M_fractures->getFracture ( f )->getData().getUl();
	        		//std::cout << " M_fractures->getFracture ( f )->getData().getUl() " << M_fractures->getFracture ( f )->getData().getUl() << std::endl;
	        		u[ n-1 ] =M_fractures->getFracture ( f )->getData().getUr();

	        }
	        else if ( M_fractures->getFracture ( f )->getData().getBc() =="p")
	        {
	        	u[ 0 ] = u[ n-2 ];
	        	u[ n-1 ] = u[ 1 ];

	        }
	        else
	        {
	        	std::cout << " bc non è stato impostato " << std::endl;
	        }

			u0 = u;
		}

		// esporto la soluzione
	//	std::ostringstream ss;
	//	ss << i;
	//	std::string name("u" );
	//	name = name + ss.str();
		std::ofstream exp("u.txt");

		exp << u;


	}

	return;
}// solve


