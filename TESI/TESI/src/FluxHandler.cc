
#include "../include/FluxHandler.h"


/**************************************************************************/
/*  FluxHandler.cc 	                                                	  */
/*  Classe che rappresenta una funzione flusso					 		  */
/**************************************************************************/

FluxHandler::FluxHandler ( const GetPot& dataFile, const std::string& sectionFlux ):
						 M_sectionFlux ( sectionFlux ),
						 M_flux ( dataFile ( ( M_sectionFlux + "f" ).data (), "x*x./2" ) ),
						 M_flux1 ( dataFile ( ( M_sectionFlux + "f1" ).data (), "x" ) ),
						 M_Monotone ( dataFile ( ( M_sectionFlux + "monotone" ).data (), "true" ) ),
						 M_H ( dataFile ( ( M_sectionFlux + "h" ).data (), 1.0 ) ),
						 M_Us ( dataFile ( ( M_sectionFlux + "us" ).data (), 0.5 ) ),
						 M_first ( dataFile ( ( M_sectionFlux + "us_meno" ).data (), 0.5 ) ),
						 M_second ( dataFile ( ( M_sectionFlux + "us_piu" ).data (), 0.5 ) ),
						 M_a ( dataFile ( ( M_sectionFlux + "bc" ).data (), 1. ) ),
						 M_b ( dataFile ( ( M_sectionFlux + "bc" ).data (), 1. ) )
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

