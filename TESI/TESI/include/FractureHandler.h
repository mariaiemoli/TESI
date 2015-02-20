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
#include "BCHandler.h"

class FractureHandler
{
public:

	/**
	 * Flags che uso per distinguere sulla mesh di ogni level set le regioni " non tagliate ",
	 * cioè dove non vi è intersezione con altre fratture, dalle regioni " tagliate ", cioè dove il level set ha un'intersezione.
	 */
    enum
    {
        FRACTURE_UNCUT = 10000,
        FRACTURE_INTERSECT = 10000
    };


    FractureHandler ( const GetPot& dataFile, const size_type ID,
    				  const ExporterPtr_Type& exporter,
    				  const std::string& section = "fractureData/" );

    void init ();

    /**
     * Funzione che calcola il vettore normale alla frattura e la mappa di conversione dalla mesh reale alla mesh piatta.
     */
    void normalVectorAndMap ( const getfem::mesh_fem& mediumMeshFEMPressure );


    sizeVector_Type getDofIntersection( );


    /**
     * Funzione che calcola il passo di griglia, h^-1.
     */
    void computeInvH ( const BCHandlerPtr_Type& bcHandler );


    inline const size_type& getID() const
	{
		return M_ID;
	}

    void numFractures ( const size_type& numFractures )
    {
        M_meshLevelSetIntersect.resize ( numFractures );
        M_levelSetIntersect.resize ( numFractures );
        M_fractureIntersectElements.resize ( numFractures );
        M_fractureIntersectElementsGlobalIndex.resize ( numFractures );
    }

    /**
     * Funzione che, data un'altra frattura con cui si interseca, imposta i  valori legati all'intersezione: costruisce sulla mesh
     * la regione " tagliata ", aggiunge i gradi di libertà estesi,  e imposta che tali gradi di libertà siano gli stessi sulle due fratture.
     */
    void setMeshLevelSetFracture ( FractureHandler& otherFracture );//, size_type& globalIndex, const std::string& type );


    inline FractureData& getData()
    {
    	return M_data;
    }

    inline LevelSetHandlerPtr_Type& getLevelSet()
    {
    	return M_levelSet;
    }

    inline scalar_type getH () const
    {
    	return M_h;
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

    /**
     * \return M_LevelSetIntersect[f]: level set di indice f intersecato dalla frattura corrente
     */
    GFLevelSetPtr_Type getLevelSetIntersect ( const size_type& f )
    {
        return M_levelSetIntersect[f];
    }

    inline const bgeot::pgeometric_trans& getGeometricTransformation ( ) const
    {
        return M_geometricTransformation;
    }

    inline const getfem::mesh_fem& getMeshFEM ( ) const
    {
        return M_meshFEM;
    }

    inline const getfem::mesh_fem& getMeshFEM2 ( ) const
    {
        return M_meshFEM2;
    }

    inline const getfem::mesh_im& getIntegrationMethodVisualization ( ) const
    {
        return M_integrationMethodVisualization;
    }

    inline const getfem::mesh_im& getIntegrationMethod ( ) const
    {
        return M_integrationMethod;
    }

    inline const getfem::mesh_im& getIntegrationMethod2 ( ) const
    {
        return M_integrationMethod2;
    }

    inline const getfem::mesh_fem& getMeshFEMVisualization ( ) const
    {
        return M_meshFEMVisualization;
    }

    inline const scalarVector_Type& getMagnificationMapFactor1 ( ) const
    {
        return M_magnificationMapFactor1;
    }


    inline const scalar_type& getMagnificationMapFactor1 ( const size_type& dof ) const
    {
        return M_magnificationMapFactor1 [ dof ];
    }


    inline const scalarVector_Type& getMagnificationMapFactor2 ( ) const
    {
        return M_magnificationMapFactor2;
    }


    inline const scalar_type& getMagnificationMapFactor2 ( const size_type& dof ) const
    {
        return M_magnificationMapFactor2 [ dof ];
    }

    inline const scalarVector_Type& getEtaTangentialInterpolated ( ) const
    {
        return M_etaTangentialInterpolated;
    }


    inline const scalar_type& getEtaTangentialInterpolated ( const size_type& dof ) const
    {
        return M_etaTangentialInterpolated [ dof ];
    }

    inline const scalarVector_Type& getInverseMeshSize ( ) const
    {
        return M_inverseMeshSize;
    }


    inline const scalar_type& getInverseMeshSize ( const size_type& dof ) const
    {
        return M_inverseMeshSize [ dof ];
    }

    const pairSizeVectorContainer_Type& getFractureIntersectElementsGlobalIndex () const
    {
        return M_fractureIntersectElementsGlobalIndex;
    }

    pairSizeVectorContainer_Type& getFractureIntersectElementsGlobalIndex ()
    {
        return M_fractureIntersectElementsGlobalIndex;
    }


    inline scalarVector_Type getCi () const
    {
    	return M_ci;
    }

    //scalarVector_Type getDofIntersection () const;


private:

    size_type M_ID;

    FractureData M_data;

    LevelSetHandlerPtr_Type M_levelSet;

	scalar_type M_h;

    ExporterPtr_Type M_exporter;

    getfem::mesh M_mesh;
    getfem::mesh M_meshFlat;

    GFMeshLevelSetPtrContainer_Type M_meshLevelSetIntersect;
    GFLevelSetPtrContainer_Type M_levelSetIntersect;

    sizeVectorContainer_Type M_fractureIntersectElements;
    pairSizeVectorContainer_Type M_fractureIntersectElementsGlobalIndex;

    // eta_gamma = d/K_normale - vector
    scalarVector_Type M_etaNormalInterpolated;
    // eta_t=1/(K_t*d) - vector
    scalarVector_Type M_etaTangentialInterpolated;

    getfem::mesh_fem M_meshFEM;
    getfem::mesh_fem M_meshFEM2;

    getfem::mesh_fem M_meshFEMVisualization;

    // integration method
    getfem::mesh_im M_integrationMethod;
    // integration method
    getfem::mesh_im M_integrationMethodVisualization;
    getfem::mesh_im M_integrationMethod2;

    //fattore di scala nella mappatura fra frattura piana e mappata
    scalarVector_Type M_magnificationMapFactor1;
    //fattore di scala nella mappatura fra frattura piana e mappata
    scalarVector_Type M_magnificationMapFactor2;

    //componenti della normale
    scalarVector_Type M_normal1;
    scalarVector_Type M_normal2;

    // M_mediumInverseMeshSize = 1.0 / M_mediumMeshSize; per la frattura
    scalarVector_Type M_inverseMeshSize;

    scalarVector_Type M_ci;

    // integration method (velocity)
    getfem::mesh_im M_integrationMethodLinear;
    // mesh_fem for coefficients
    getfem::mesh_fem M_meshFEMLinear;


    // Geometric transformation usign pressure finite elements type
    bgeot::pgeometric_trans M_geometricTransformation;

};


typedef FractureHandler FractureHandler_Type;											/*!< classe FractureHandler */
typedef boost::shared_ptr<FractureHandler_Type> FractureHandlerPtr_Type;				/*!< puntatore alla classe FractureHandler */
typedef std::vector<FractureHandlerPtr_Type> FracturePtrContainer_Type;					/*!< vettore di puntatori alla classe FractureHandler */
typedef std::vector<FractureHandler_Type> FractureContainer_Type;						/*!< vettore di classi FractureHandler */
typedef boost::shared_ptr<FracturePtrContainer_Type> FracturePtrContainerPtr_Type;		/*!< puntatore a un vettore di puntatori alla classe FractureHandler */


#endif /* _FRACTUREHANDLER_H_ */
