
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
						 M_Uin ( dataFile ( ( M_sectionFlux + "uin" ).data (), 0. ) ),
						 M_Uout ( dataFile ( ( M_sectionFlux + "uout" ).data (), 0. ) )
{}


void FluxHandler::update_Bc ( const size_type& pos, const scalar_type& u )
{
	if ( pos == 0 )
	{
		M_Uin = u;
	}
	else
	{
		M_Uout = u;
	}

	return;
}// update_Bc
