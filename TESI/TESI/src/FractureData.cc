
#include "../include/FractureData.h"

/**************************************************************************/
/*  FractureData.cc                                              	   	  */
/*  Classe che contiene i dati per il problema di saturazione	 		  */
/**************************************************************************/

FractureData::FractureData ( const GetPot& dataFile,
				   	   	   	 const std::string& section,
			                 const std::string& sectionDomain,
				   	   	   	 const std::string& sectionSaturation ):
				   	   	   	 M_section ( section ),
				   	   	   	 M_sectionDomain ( M_section + sectionDomain ),
				   	   	   	 M_sectionSaturation ( M_section + sectionSaturation ),
				   	   	   	 // domain
				             M_thickness ( dataFile ( ( M_sectionDomain + "thickness" ).data (), 0.01 ) ),
				             M_lengthAbscissa ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
				             M_lengthOrdinate ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 0. ) ),
				             M_lengthQuota ( dataFile ( ( M_sectionDomain + "lengthQuota" ).data (), 0. ) ),
				             M_a ( dataFile ( ( M_sectionDomain + "a" ).data (), -1. ) ),
							 M_b ( dataFile ( ( M_sectionDomain + "b" ).data (), 1. ) ),
				             M_spatialDiscretization ( dataFile ( ( M_sectionDomain + "spatialDiscretization" ).data (), 200 ) ),
				             M_translateAbscissa ( dataFile ( ( M_sectionDomain + "translateAbscissa" ).data (), 0. ) ),
				             M_spaceDimension ( dataFile ( ( M_section + "spaceDimension" ).data (), 1. ) ),
				             M_meshType ( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(1,1)" ) ),
				             M_meshSpacing( dataFile ( ( M_sectionDomain + "spacing" ).data (), "x" ) ),
				             M_fEMType ( dataFile ( ( M_sectionDomain + "FEMType" ).data (), "FEM_PK(1,0)" ) ),
				             M_fEMType2 ( dataFile ( ( M_sectionDomain + "FEMType2" ).data (), "FEM_PK(1,1)" ) ),
				             M_integrationType ( dataFile ( ( M_sectionDomain + "integrationType" ).data (),
				                                                    "IM_GAUSS1D(2)" ) ),
							 M_integrationType2 ( dataFile ( ( M_sectionDomain + "integrationType2" ).data (),
																    "IM_GAUSS1D(3)" ) ),
				             // saturation
				   	   	   	 M_u0 ( dataFile ( ( M_sectionSaturation + "u0").data (), "1." ) ),
				   	   	   	 M_numberFlux( dataFile ( ( M_sectionSaturation + "numberFlux").data (), 1 ) ),
				   	   	   	 M_x( dataFile ( ( M_sectionSaturation + "x0").data (), 0. ) ),
				             M_invP ( dataFile ( (M_sectionSaturation + "invPorosity" ).data (), "1." ) )
{
	M_flux.resize( M_numberFlux );

	std::ostringstream sectionFlux;

	for ( size_type i = 0; i < M_numberFlux; i++ )
	{
		sectionFlux << M_sectionSaturation << "fluxData" << i << "/";

		M_flux [ i ].reset( new FluxHandler_Type( dataFile, sectionFlux.str() ));

		sectionFlux.str("");
	}

	M_Int = 1000000;

}// costruttore


void FractureData::feval( scalarVector_Type& flux, const scalarVector_Type& u0, const size_type i, const size_type f )
{
    flux.clear();

    std::string Flux;

    for ( size_type j = 0; j < u0.size(); j++ )
    {
		scalar_type P = porosity ( u0 [ j ]);

		if( i == 1)
		{
			Flux = M_flux [ f ]->getFlux();
		}
		else
		{
			Flux = M_flux [ f ]->getFlux1();
		}

		M_parser.setString ( Flux );
    	M_parser.setVariable ( "x", u0 [ j ] );

		flux.push_back( P*M_parser.evaluate () );

    }

    return;
}// feval


scalar_type FractureData::feval_scal( const scalar_type& us, const size_type i, const size_type f )
{
	//	scalar_type P = porosity ( us );

	std::string Flux;

	if ( f == 0 )
	{
		Flux = M_flux [ i ]->getFlux();

		M_parser.setString ( Flux );

		M_parser.setVariable ( "x", us );

		//return P*M_parser.evaluate ();

		return M_parser.evaluate ();
	}
	else
	{
		Flux = M_flux [ i ]->getFlux1();

		M_parser.setString ( Flux );

		M_parser.setVariable ( "x", us );

		//return P*M_parser.evaluate ();

		return M_parser.evaluate ();


	}

}// feval_scal


std::string FractureData::getFlux( const size_type& i )
{
	return M_flux [ i ]->getFlux();
}// getFlux


std::string FractureData::getFlux1( const size_type& i )
{
	return M_flux [ i ]->getFlux1();
}// getFlux1


scalar_type FractureData::meshSpacing( const scalar_type& x )
{
    M_parser.setString ( M_meshSpacing );
    M_parser.setVariable ( "x", x );

    return M_parser.evaluate ();
}// meshSpacing


scalar_type FractureData::porosity( const scalar_type& u )
{
	M_parser.setString( M_invP );
	M_parser.setVariable ( "x", u );

	return M_parser.evaluate ();
}// porosity


scalar_type FractureData::getCi ( const scalar_type& x )
{

	M_parser.setString( M_u0 );
	M_parser.setVariable ( "x", x );

	return M_parser.evaluate ();
}// getCi



void FractureData::update_UI ( const size_type& j, const scalar_type& u )
{

	for ( size_type i = 0; i < M_flux.size(); i++ )
	{
		M_flux [ i ]-> update_UI ( j, u );
	}

	return;
}// update_UI


void FractureData::updateSI ( const scalarVector_Type& f )
{
	if ( M_Int != 1000000 )
	{
		M_SI = f [ M_Int ];
	}
	else
	{
		M_SI = f [ M_spatialDiscretization ];
	}

	return;
}// updateSI


void FractureData::imposeIntersection ( const size_type& i )
{
	M_Int = i;
}// imposeIntersection
