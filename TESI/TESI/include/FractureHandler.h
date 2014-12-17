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

#ifndef _FRACTUREHANDLER_H_
#define _FRACTUREHANDLER_H_ 1

/**************************************************************************/
/*  FractureHandler.h                                                 	  */
/*  Classe che inizializza e gestisce una frattura				 		  */
/**************************************************************************/

#include "StringUtility.h"
#include "Core.h"
#include "FractureData.h"
#include "LevelSetHandler.h"
#include "UsefulFunctions.h"
#include "Exporter.h"

class FractureHandler
{
public:

    FractureHandler ( const GetPot& dataFile, const size_type ID,
    				  const ExporterPtr_Type& exporter,
    				  const std::string& section = "fractureData/",
    				  const std::string& sectionSaturation = "saturation/" );

    void init ();

    /**
     * Funzione che definisce il nodo in cui la frattura corrente ha un'intersezione, e la frattura coinvolta
     */
    void setFractureIntersection ( const sizeVector_Type& nodes, const FractureHandler& fractureInvolved );


    inline const size_type& getID() const
	{
		return M_ID;
	}

    inline FractureData& getData()
    {
    	return M_data;
    }

    inline LevelSetHandlerPtr_Type& getLevelSet()
    {
    	return M_levelSet;
    }

    inline const getfem::mesh& getMesh ( ) const
    {
        return M_mesh;
    }


    inline getfem::mesh& getMesh ( )
    {
        return M_mesh;
    }

    inline const getfem::mesh& getMeshFlat ( ) const
    {
        return M_meshFlat;
    }


    inline getfem::mesh& getMeshFlat ( )
    {
        return M_meshFlat;
    }

    inline scalarVector_Type getCi () const
    {
    	return M_ci;
    }

    void numFractures ( const size_type& numFractures )
    {
        M_fractureIntersectElements.resize ( numFractures );
        M_fractureIntersectElementsGlobalIndex.resize ( numFractures );
    }



private:

    size_type M_ID;

    FractureData M_data;

    LevelSetHandlerPtr_Type M_levelSet;

    ExporterPtr_Type M_exporter;

    getfem::mesh M_mesh;
    getfem::mesh M_meshFlat;

    scalarVector_Type M_ci;

    // Geometric transformation usign pressure finite elements type
    bgeot::pgeometric_trans M_geometricTransformation;

    sizeVectorContainer_Type M_fractureIntersectElements;
    pairSizeVectorContainer_Type M_fractureIntersectElementsGlobalIndex;  // vettore di vettori di coppie

};


typedef FractureHandler FractureHandler_Type;											/*!< classe FractureHandler */
typedef boost::shared_ptr<FractureHandler_Type> FractureHandlerPtr_Type;				/*!< puntatore alla classe FractureHandler */
typedef std::vector<FractureHandlerPtr_Type> FracturePtrContainer_Type;					/*!< vettore di puntatori alla classe FractureHandler */
typedef std::vector<FractureHandler_Type> FractureContainer_Type;						/*!< vettore di classi FractureHandler */
typedef boost::shared_ptr<FracturePtrContainer_Type> FracturePtrContainerPtr_Type;		/*!< puntatore a un vettore di puntatori alla classe FractureHandler */


#endif /* _FRACTUREHANDLER_H_ */
