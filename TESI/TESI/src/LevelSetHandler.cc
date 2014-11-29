
# include "../include/LevelSetHandler.h"

/**************************************************************************/
/*  LevelSetHandler.cc                                                	  */
/*  Classe che inizializza e gestisce il level set				 		  */
/**************************************************************************/

LevelSetHandler::LevelSetHandler ( const GetPot& dataFile,
								   const std::string& section,
								   const std::string& sectionLevelSet ):
								   M_Data ( new LevelSetData_Type ( dataFile, section, sectionLevelSet ) )
{}

