
# include "../include/BC.h"

/**************************************************************************/
/*  BC.cc															      */
/*  Classe che introduce le condizioni al bordo sul problema              */
/**************************************************************************/


BC::BC ( const GetPot& dataFile,
		 const std::string& section,
		 getfem::mesh& mesh,
		 const std::string& MeshType,
		 const std::string& sectionSaturation ):
		 M_meshFEM(mesh),
		 M_section ( section ),
		 M_sectionSaturation ( M_section + sectionSaturation )
{};
