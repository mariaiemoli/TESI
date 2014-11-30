/*
 * Progetto di tesi
 * Politecnico di Milano
 *
 * \author Iemoli Maria
 * \date 2015
 *
 * Problema di saturazione per un fluido bifase in un network di fratture
 *
 */

#include "include/SaturationFractured.h"
#include "include/FracturesSet.h"
#include "include/Core.h"

/**************************************************************************/
/*  main program                                                          */
/**************************************************************************/


int main ( int argc, char* argv [ ] )
{
	std::string fileName("data");

	if ( argc == 2 )
	{
		fileName = argv[1];
	}

	GetPot dataFile( fileName.c_str() );

	const std::string section = "";

	const std::string vtkFolder = "vtk/";


	//Data exporter
	std::cout << "Create the data exporter..." << std::flush;
	ExporterPtr_Type exporter( new Exporter_Type(dataFile));
	std::cout << " completed!" <<std::endl;


	// Creo la cartella dove salvare i risultati se già non esiste
//	std::string s = "mkdir "+vtkFolder;
//	system (s.c_str());


	// Fracture Set
	std::cout << "Create the set of fractures for " << std::flush;
	const size_type numberFractures = dataFile(
			(section + "numberFractures").data(), 0);
	std::cout << numberFractures << " fracture(s)..." << std::endl << std::flush;


	// Inizializzo le fratture
	FracturesSetPtr_Type fractures ( new  FracturesSet  );

	std::cout << " completed! " << std::endl;
	fractures->init( dataFile, section, numberFractures, exporter );


	// Risolvo il problema in saturazione
	SaturationFracturedPtr_Type saturation(new SaturationFractured_Type( dataFile, fractures ) );

	std::cout << " risolvo il problema " << std::endl;

	saturation->init();

	saturation-> solve();

	return 0;

}
