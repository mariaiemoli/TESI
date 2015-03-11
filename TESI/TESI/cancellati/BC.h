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


#ifndef BC_H_
#define BC_H_ 1

#include "Core.h"

/**************************************************************************/
/*  BC.h															      */
/*  Classe che introduce le condizioni al bordo sul problema              */
/**************************************************************************/

class BC
{
public:

	// Costruttore
    BC ( const GetPot& dataFile,
    	 const std::string& section,
    	 getfem::mesh& mesh,
         const std::string& MeshType,
		 const std::string& sectionSaturation = "saturation/");

	/*
    inline const sizeVector_Type& getDirichlet ( ) const
    {
        return M_dirichlet;
    }


    inline const size_type& getDirichlet ( const size_type& dof ) const
    {
        return M_dirichlet [ dof ];
    }


    inline const sizeVector_Type& getNeumann ( ) const
    {
        return M_neumann;
    }


    inline const size_type& getNeumann ( const size_type& dof ) const
    {
        return M_neumann [ dof ];
    }


    inline const getfem::mesh_fem& getMeshFEM ( ) const
    {
        return M_meshFEM;
    }
*/

private:


    getfem::mesh_fem M_meshFEM;

	std::string M_section;
	std::string M_sectionSaturation;


    // flags for BC
    sizeVector_Type M_dirichlet;
    sizeVector_Type M_neumann;
    sizeVector_Type M_extBoundary;

};

typedef BC BC_Type;												/*!< Classe BC */
typedef boost::shared_ptr<BC> BCPtr_Type;						/*!< Puntatore alla classe BC */
typedef std::vector<BCPtr_Type> BCPtrContainer_Type;			/*!< Vettore di puntatori alla classe BC */

#endif /* BC_H_ */
