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

class FractureData
{
public:

	FractureData ( const GetPot& dataFile,
				   const std::string& section = "fractureData/",
                   const std::string& sectionDomain = "domain/",
				   const std::string& sectionSaturation = "saturation/" );

	void feval( scalarVector_Type& fl, const scalarVector_Type& u0, const size_type i );

	scalar_type feval_scal( const scalar_type& us);

	scalar_type fzero( const std::string& f, const scalar_type& a, const scalar_type& b );

    scalar_type meshSpacing ( const scalar_type& x );

    scalar_type velocity( const scalar_type& u );

    scalar_type porosity( const scalar_type& u );


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

    inline scalar_type getH () const
    {
    	return M_h;
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

	inline std::string getBc () const
	{
		return M_bc;
	}

	inline scalar_type getUl () const
	{
		return M_ul;
	}

	inline scalar_type getUr () const
	{
		return M_ur;
	}

	inline scalar_type getX0 () const
	{
		return M_x0;
	}

	inline std::string getFlux()
	{
		return M_flux;
	}

	inline std::string getFlux1()
	{
		return M_flux1;
	}

	inline scalar_type getCfl()
	{
		return M_cfl;
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

    size_type M_spatialDiscretization;
    scalar_type M_translateAbscissa;
    bgeot::dim_type M_spaceDimension;

	std::string M_meshType;
    std::string M_meshSpacing;

    // saturation
	std::string M_bc;

	scalar_type M_ul;
	scalar_type M_ur;

	scalar_type M_x0;

	std::string M_flux;
	std::string M_flux1;

	scalar_type M_cfl;

	std::string M_U;
	std::string M_invP;

	mutable LifeV::Parser M_parser;

};


#endif /* _FRACTUREDATA_H_ */
