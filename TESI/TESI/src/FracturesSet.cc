
#include "../include/FracturesSet.h"

/**************************************************************************/
/*  FracturesSet.cc                                                 	  */
/*  Classe che rappresenta l'insieme di tutte le fratture		 		  */
/**************************************************************************/

FracturesSet::FracturesSet():
				M_intersections(new FractureIntersect_Type)
{}// costruttore nullo


void FracturesSet::init ( const GetPot& dataFile, const std::string& section,
						  const size_type& numFractures, getfem::mesh& mesh, getfem::mesh_level_set& meshLevelSet,
						  const std::string& integrationTypeVelocity, const getfem::mesh_fem& meshFEMScalar,
						  const getfem::mesh_fem& meshFEMVector,
						  const ExporterPtr_Type& exporter )
{
	M_fractures.resize( numFractures );

	std::ostringstream sectionFracture;

	for ( size_type f = 0; f < numFractures; f++ )
	{
		sectionFracture << section << "fractureData" << f << "/";

		M_fractures [ f ].reset ( new FractureHandler_Type ( dataFile, f, exporter, sectionFracture.str() ));

		M_fractures [ f ]->init();

		M_fractures [ f ]->numFractures ( numFractures );

		// Per ogni frattura inizializzo il level set
		M_fractures [ f ]->getLevelSet()->init ( mesh, integrationTypeVelocity, meshFEMScalar, meshFEMVector );

		// Calcolo la normale alla frattura e il fattore di conversione della lunghezza della mesh reale
		M_fractures [ f ]->normalVectorAndMap ( meshFEMScalar );

		// Per ogni frattura aggiungo le informazioni relative al corrispondente level set a meshLevelSet
		meshLevelSet.add_level_set ( M_fractures[f]->getLevelSet()->getLevelSet() );

		sectionFracture.str("");
	}

	// Devo aggiornare la mesh, sono stati aggiunti i levelset
	meshLevelSet.adapt ();

	// costruisco l'intersezione
	M_intersections->constructIntesection ( M_fractures );

	return;
}// init
