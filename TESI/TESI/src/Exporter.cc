
#include "../include/Exporter.h"

/**************************************************************************/
/*  Exporter.h                                                 			  */
/*  Libreria in cui definisco le funzioni per esportare i dati			  */
/**************************************************************************/

Exporter::Exporter ( const GetPot& dataFile, const std::string& section ) :
					M_vtkFolder(dataFile((section + "folderVTK").data(), "./vtk/"))
{
}// costruttore


void Exporter::spy ( const sparseMatrixPtr_Type& matrix, const std::string& nameFile ) const
{
    gmm::MatrixMarket_IO::write(nameFile.c_str(), *matrix);

    return;

} // spy


void Exporter::spy ( const scalarVectorPtr_Type& vector, const std::string& nameFile ) const
{
    std::ofstream file;
    file.open(nameFile.c_str());
    for ( size_type i = 0; i < vector->size(); ++i )
    {
        file << std::setprecision(15) << vector->at(i) << std::endl;
    }
    file.close();

    return;

} // spy


