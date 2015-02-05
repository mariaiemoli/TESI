/* Progetto di tesi
 * Politecnico di Milano
 *
 * \author Iemoli Maria
 * \date 2015
 *
 * Problema di saturazione per un fluido bifase in un network di fratture
 *
 */

#ifndef _FRACTUREDATA_H_
#define _FRACTUREDATA_H_ 1

/**************************************************************************/
/*  FractureData.h                                              	   	  */
/*  Classe che contiene i dati per il problema di saturazione	 		  */
/**************************************************************************/

#include "StringUtility.h"
#include "Core.h"
#include "Parser.h"
#include "FluxHandler.h"

class FractureData
{
public:

	FractureData ( const GetPot& dataFile,
				   const std::string& section = "fractureData/",
                   const std::string& sectionDomain = "domain/",
				   const std::string& sectionSaturation = "saturation/" );

	void feval( scalarVector_Type& fl, const scalarVector_Type& u0, const size_type i, const size_type f );

	scalar_type feval_scal( const scalar_type& us, const size_type i, const size_type f = 0 );

	std::string getFlux( const size_type& i );

	std::string getFlux1( const size_type& i );

    scalar_type meshSpacing ( const scalar_type& x );

    scalar_type porosity( const scalar_type& u );

    scalar_type getCi ( const scalar_type& x );

    /**
     * Funzione che aggiorna il valore della saturazione nel nodo di intersezione
     */
    void update_UI ( const size_type& i, const scalar_type& u );

    /**
     * Funzione che agggiorna il valore della funzione flusso all'intersezione
     */
    void updateSI ( const scalarVector_Type& f );

    void imposeIntersection ( const size_type& i );

    inline scalar_type getThickness () const
    {
        return M_thickness;
    }

	inline scalar_type getLengthAbscissa () const
	{
		return M_lengthAbscissa;
	}

	inline scalar_type getLengthOrdinate () const
	{
		return M_lengthOrdinate;
	}

	inline scalar_type getLengthQuota () const
	{
		return M_lengthQuota;
	}

    inline scalar_type getA () const
    {
    	return M_a;
    }

    inline scalar_type getB () const
    {
    	return M_b;
    }

	inline scalar_type getSpatialDiscretization () const
	{
		return M_spatialDiscretization;
	}

	inline scalar_type getSpaceDimension () const
	{
		return M_spaceDimension;
	}

    inline std::string getMeshType () const
    {
        return M_meshType;
    }

    inline std::string getFEMType () const
    {
        return M_fEMType;
    }

    inline std::string getFEMType2 () const
    {
        return M_fEMType2;
    }

    inline std::string getMeshSpacing () const
    {
        return M_meshSpacing;
    }


    inline std::string getIntegrationType () const
    {
        return M_integrationType;
    }

    inline std::string getIntegrationType2 () const
    {
        return M_integrationType2;
    }

	inline size_type getNumberFlux () const
	{
		return M_numberFlux;
	}

	inline FluxHandlerPtr_Type getFluxHandler( const size_type i )
	{
		return M_flux[ i ];
	}

	inline scalar_type getXd () const
	{
		return M_x;
	}

	inline scalar_type getSI()
	{
		return M_SI;
	}

private:

	std::string M_section;
    std::string M_sectionDomain;
	std::string M_sectionSaturation;

	// domain
    scalar_type M_thickness; // thickness of the fracture

    scalar_type M_lengthAbscissa;
	scalar_type M_lengthOrdinate;
	scalar_type M_lengthQuota;
	scalar_type M_a;
	scalar_type M_b;
	scalar_type M_h;

	scalar_type M_spatialDiscretization;

    scalar_type M_translateAbscissa;
    bgeot::dim_type M_spaceDimension;

	std::string M_meshType;
    std::string M_meshSpacing;

    std::string M_fEMType;
    std::string M_fEMType2;

    std::string M_integrationType;
    std::string M_integrationType2;

    // saturation
	std::string M_u0;

	size_type M_numberFlux;

	// punto in cui c'è l'eventuale discontinuità nella funzione flusso
	scalar_type M_x;

	FluxPtrContainer_Type M_flux;
	size_type M_Int;
	scalar_type M_SI;

	std::string M_invP;

	mutable LifeV::Parser M_parser;

};


#endif /* _FRACTUREDATA_H_ */
