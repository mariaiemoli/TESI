
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
				             // saturation
				   	   	   	 M_bc ( dataFile ( ( M_sectionSaturation + "bc").data (), "p" ) ),
				   	   	   	 M_u0 ( dataFile ( ( M_sectionSaturation + "u0").data (), "1." ) ),
				   	   	   	 M_numberFlux( dataFile ( ( M_sectionSaturation + "numberFlux").data (), 1 ) ),
				   	   	   	 M_x( dataFile ( ( M_sectionSaturation + "x0").data (), 0. ) ),
				             M_cfl ( dataFile ( ( M_sectionSaturation + "cfl" ).data (), 1.8 ) ),
				             M_U ( dataFile ( (M_sectionSaturation + "U" ).data (), "1." ) ),
				             M_invP ( dataFile ( (M_sectionSaturation + "invPorosity" ).data (), "1." ) )
{
	M_h = ( M_b-M_a )/(M_spatialDiscretization-2);

	M_flux.resize( M_numberFlux );

	std::ostringstream sectionFlux;

	for ( size_type i = 0; i < M_numberFlux; i++ )
	{
		sectionFlux << M_sectionSaturation << "fluxData" << i << "/";

		M_flux [ i ].reset( new FluxHandler_Type( dataFile, sectionFlux.str() ));

	}

}// costruttore


void FractureData::feval( scalarVector_Type& flux, const scalarVector_Type& u0, const size_type i, const size_type f )
{
    flux.clear();

    std::string Flux;

    for ( size_type j = 0; j < u0.size(); j++ )
    {
		scalar_type U = velocity ( u0 [ j ]);
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

        flux.push_back( U*P*M_parser.evaluate () );

    }
    return;
}// feval


// farla con un template
scalar_type FractureData::feval_scal( const scalar_type& us, const size_type f )
{
	scalar_type U = velocity ( us );
	scalar_type P = porosity ( us );

	std::string Flux;

	if ( f == 0 )
	{
		Flux = M_flux [ 0 ]->getFlux();

		M_parser.setString ( Flux );

		M_parser.setVariable ( "x", us );

		return U*P*M_parser.evaluate ();
	}
	else
	{
		Flux = M_flux [ 0 ]->getFlux1();

		M_parser.setString ( Flux );

		M_parser.setVariable ( "x", us );

		return U*P*M_parser.evaluate ();

	}

}// feval_scal

scalar_type FractureData::fzero( const std::string& f, const scalar_type& a, const scalar_type& b )
{
	scalar_type toll = 1.0E-7;
	scalar_type x0;

	M_parser.setString ( f );

    M_parser.setVariable ( "x", a );

	scalar_type f_a = M_parser.evaluate ();

	if ( f_a < toll )
	{
		return a;
	}

	M_parser.setVariable ( "x", b );

	scalar_type f_b = M_parser.evaluate ();

	if ( f_b < toll )
	{
		return b;
	}

	scalar_type c = ( a + b )/2;

	scalar_type f_c = M_parser.evaluate ();

	if ( f_b < toll )
	{
		return b;
	}

	if ( f_a * f_c <0 )
	{
		x0 = fzero( f, a, c );
	}
	else
	{
		x0 = fzero( f, c, b );
	}

	return x0;
}// fzero


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


scalar_type FractureData::velocity( const scalar_type& u )
{
	M_parser.setString( M_U );
	M_parser.setVariable ( "x", u );

	return M_parser.evaluate ();
}// velocity


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

	/*
	if( x == -4.0 )
	{
		return 1.;
	}
	else
	{
		return 0;
	}
	*/
}// getCi
