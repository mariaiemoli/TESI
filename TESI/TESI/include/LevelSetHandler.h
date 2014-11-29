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

#ifndef _LEVELSETHANDLER_H_
#define _LEVELSETHANDLER_H_ 1

#include "LevelSetData.h"

/**************************************************************************/
/*  LevelSetHandler.h                                                 	  */
/*  Classe che inizializza e gestisce il level set				 		  */
/**************************************************************************/

class LevelSetHandler
{
public:
	LevelSetHandler ( const GetPot& dataFile, const std::string& section =  "fractureData/",
			  	  	  const std::string& sectionLevelSet = "levelSet/" );


    inline LevelSetDataPtr_Type& getData ( )
    {
        return M_Data;
    }

private:

	LevelSetDataPtr_Type M_Data;

    // M_mediumMesh, level set
    GFMeshLevelSetPtr_Type M_mesh;
    // Il level set che definisce la frattura
    GFLevelSetPtr_Type M_levelSet;


};


typedef LevelSetHandler LevelSetHandler_Type;									/*!< Classe LevelSetHandler */
typedef boost::shared_ptr<LevelSetHandler_Type> LevelSetHandlerPtr_Type;		/*!< Puntatore alla classe LeverSetHandler */


#endif /* _LEVELSETHANDLER_H_ */
