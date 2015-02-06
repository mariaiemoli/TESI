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
#include "include/MediumData.h"
//#include "include/DarcyFractured.h"
#include "include/MeshHandler.h"



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

	// Creo la cartella dove salvare i risultati se già non esiste
	//	std::string s = "mkdir "+vtkFolder;
	//	system (s.c_str());


	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;

	std::cout << "Solving the Saturation problem for a set of fractures." << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	//Data exporter
	std::cout << "Create the data exporter..." << std::flush;
	ExporterPtr_Type exporter( new Exporter_Type( dataFile ));
	std::cout << " completed!" <<std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	//Mesh Handler
	std::cout << "Create the meshHandler..." << std::flush;
	MeshHandlerPtr_Type mesh(new MeshHandler_Type( dataFile, "mediumData/domain/" ));
	mesh->setUpMesh();
	mesh->setUpFEM();
	std::cout<< " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;



	// Medium data for the Darcy problem
	std::cout << "Create the mediumData for the Darcy problem.." << std::flush;
	const std::string sectionSolverDarcy = "darcy/";
	MediumDataPtr_Type mediumDataDarcy( new MediumData_Type( dataFile, sectionSolverDarcy ));
	std::cout << " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	// Fracture Set
	std::cout << "Create the set of fractures for " << std::flush;
	const size_type numberFractures = dataFile((section + "numberFractures").data(), 0 );
	std::cout << numberFractures << " fracture(s)..."  << std::flush;


	// Inizializzo le fratture
	FracturesSetPtr_Type fractures ( new  FracturesSet  );

	std::cout << " completed! " << std::endl;
	fractures->init( dataFile, section, numberFractures, exporter );
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;

/*
	// Fracture boundary conditions
	std::cout << "Create fracture boundary conditions..." << std::flush;
	BCPtrContainer_Type bcFracture(numberFractures);
	for ( size_type f = 0; f < numberFractures; ++f )
	{
		bcFracture [ f ].reset(new BC_Type( fractures->getFracture( f )->getMeshFlat(),
										   	fractures->getFracture ( f )->getData().getMeshType(),
										   	fractures->getFracture ( f ) -> getDofIntersection() ));
	}
	std::cout << " completed!" << std::endl;


	// Boundary conditions handler
	std::cout << "Create boundary conditions handler..." << std::flush;
	BCHandlerPtr_Type bcHandler(new BCHandler_Type( bcFracture ));
	std::cout << " completed!" << std::endl;
*/

	// Risolvo il problema in saturazione
	SaturationFracturedPtr_Type saturation(new SaturationFractured_Type( dataFile, fractures, exporter ) );

	saturation->init();

	saturation-> solve();

	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	return 0;

}
