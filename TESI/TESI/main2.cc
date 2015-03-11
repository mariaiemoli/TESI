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
#include "include/DarcyFractured.h"
#include "include/MeshHandler.h"



/**************************************************************************/
/*  main program                                                          */
/*  Accoppiamo completamente i problemi									  */
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

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;
	std::cout << " 			---  FLUSSO DI FLUIDI BIFASE IN MEZZI POROSI  ---			" << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;

	std::cout << " 			---  MESH DI SUPPORTO E MEZZO  ---			" << std::endl;
	std::cout << std::endl;

	//Data exporter
	std::cout << "Create the data exporter..." << std::flush;
	ExporterPtr_Type exporter( new Exporter_Type( dataFile ));
	std::cout << " completed!" <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;


	//Mesh Handler
	std::cout << "Create the meshHandler..." << std::flush;
	MeshHandlerPtr_Type mesh(new MeshHandler_Type( dataFile, "mediumData/domain/" ));
	mesh->setUpMesh();
	mesh->setUpFEM();
	std::cout<< " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;


	// Medium data for the Darcy problem
	std::cout << "Create the mediumData for the Darcy problem.." << std::flush;
	const std::string sectionSolverDarcy = "darcy/";
	MediumDataPtr_Type mediumDataDarcy( new MediumData_Type( dataFile, sectionSolverDarcy ));
	std::cout << " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	std::cout << " 			---  DEFINIZIONE FRATTURE E CONDIZIONI AL BORDO  ---			" << std::endl;
	std::cout << std::endl;

	// Fracture Set
	std::cout << "Create the set of fractures for " << std::flush;
	const size_type numberFractures = dataFile((section + "numberFractures").data(), 0 );
	std::cout << numberFractures << " fracture(s)..."  << std::flush;


	// Inizializzo le fratture
	FracturesSetPtr_Type fractures ( new  FracturesSet  );

	std::cout << " completed! " << std::endl;
	fractures->init ( dataFile, section, numberFractures, mesh->getMesh(), mesh->getMeshLevelSet(),
					  mesh->getIntegrationTypeVelocity(),
					  mesh->getMeshFEMScalar(), mesh->getMeshFEMVector(), exporter );
	std::cout << std::endl;
	std::cout << std::endl;

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
	std::cout << std::endl;
	std::cout << std::endl;


	// Boundary conditions handler
	std::cout << "Create boundary conditions handler..." << std::flush;
	BCHandlerPtr_Type bcHandler(new BCHandler_Type( bcFracture ));
	std::cout << " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;


	// Compute inverse of mesh size (h^(-1) dove h è il passo di griglia)
	std::cout << "Compute inverse of mesh size..." << std::flush;
	mesh->computeMeshMeasures();
	for ( size_type f = 0; f < numberFractures; ++f )
	{
		fractures->getFracture( f )->computeInvH(bcHandler);
	}
	std::cout << " completed!" << std::endl;
	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;



	std::cout << " 			---  PROBLEMA DI DARCY E DI SATURAZIONE ACCOPPIATI  ---			" << std::endl;
	std::cout << std::endl;

	size_type nt = 10;

	// accoppiamo i problemi
	for ( size_type i = 0; i < nt; i++ )
	{
		// Darcy problem
		DarcyFracturedPtr_Type darcy(new DarcyFractured_Type(mediumDataDarcy, mesh,
				bcHandler, fractures, exporter));

		// Initialize the solver
		darcy->init();

		// Assembly the matrices and vectors
		darcy->assembly( dataFile);

		// Solve and save the solutions
		darcy->solve();


		// Risolvo il problema in saturazione
		SaturationFracturedPtr_Type saturation(new SaturationFractured_Type( dataFile, fractures, exporter ) );

		saturation->init();

		saturation-> solve();

		// aggiorno il valore della mobilità totale per risolvere darcy
		for ( size_type f = 0; f < numberFractures; f++ )
		{
			// nel caso lineare k_{rw} = S_w
			fractures->getFracture( f )->updateEtaTangentialInterpolated( saturation-> getSaturation ( f ) );
		}
	}


	std::cout << std::endl;
	std::cout << "*******************   " << std::endl;
	std::cout << std::endl;


	return 0;

}
