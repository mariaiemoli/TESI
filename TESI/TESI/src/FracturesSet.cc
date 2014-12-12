
#include "../include/FracturesSet.h"

/**************************************************************************/
/*  FracturesSet.cc                                                 	  */
/*  Classe che rappresenta l'insieme di tutte le fratture		 		  */
/**************************************************************************/

FracturesSet::FracturesSet()//:M_intersections(new FractureIntersect_Type)
{}// costruttore nullo

void FracturesSet::init( const GetPot& dataFile, const std::string& section, const size_type& numFractures, const ExporterPtr_Type& exporter )
{
	M_fractures.resize( numFractures );

	std::ostringstream sectionFracture;

	for ( size_type f = 0; f < numFractures; f++ )
	{
		sectionFracture << section << "fractureData" << f << "/";

		M_fractures [ f ].reset ( new FractureHandler_Type ( dataFile, f, exporter, sectionFracture.str() ));

		M_fractures [ f ]->init();

		sectionFracture.str("");
	}

	// costruisco l'intersezione

	return;
}// init
