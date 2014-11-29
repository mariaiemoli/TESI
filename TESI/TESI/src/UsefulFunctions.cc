
#include "../include/UsefulFunctions.h"

/**************************************************************************/
/*  UsefulFunctions.cc                                           		  */
/*  Libreria in cui definisco funzioni ausiliarie						  */
/**************************************************************************/

void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_point_data(meshFEM, solution, solutionName.data());

    return;

}// exportSolution


void exportMesh ( const std::string& fileName, const getfem::mesh& mesh )
{
    getfem::vtk_export exporter(fileName.data(), true);

    // con true mi permette di esportare la mesh in formato txt

    exporter.exporting(mesh);

    exporter.write_mesh();

    return;
}// exportMesh
