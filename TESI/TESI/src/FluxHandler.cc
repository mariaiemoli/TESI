
#include "../include/FluxHandler.h"


/**************************************************************************/
/*  FluxHandler.cc 	                                                	  */
/*  Classe che rappresenta una funzione flusso					 		  */
/**************************************************************************/

FluxHandler::FluxHandler ( const GetPot& dataFile, const std::string& sectionFlux ):
						 M_sectionFlux ( sectionFlux ),
						 M_U ( dataFile ( ( M_sectionFlux + "U" ).data (), "v" ) ),
						 M_flux ( dataFile ( ( M_sectionFlux + "f" ).data (), "x*x./2" ) ),
						 M_flux1 ( dataFile ( ( M_sectionFlux + "f1" ).data (), "x" ) ),
						 M_Us ( dataFile ( ( M_sectionFlux + "us" ).data (), 0.5 ) ),
						 M_first ( dataFile ( ( M_sectionFlux + "us_meno" ).data (), 0.5 ) ),
						 M_second ( dataFile ( ( M_sectionFlux + "us_piu" ).data (), 0.5 ) ),
						 M_a ( dataFile ( ( M_sectionFlux + "bc1" ).data (), 1. ) ),
						 M_b ( dataFile ( ( M_sectionFlux + "bc2" ).data (), 1. ) )
{}


void FluxHandler::update_UI ( const size_type& i, const scalar_type& u )
{
	// a seconda di dove sta l'intersezione aggiorno M_a o M_b
	if ( i == 0 )
	{
		M_a = u;
	}
	else
	{
		M_b = u;
	}

	return;
}// update_Bc


void FluxHandler::monotone ( const std::string& mono )
{
	M_Monotone = mono;

	return;
}

void FluxHandler::H ( const size_type& hyp )
{
	M_H = hyp;

	return;
}

