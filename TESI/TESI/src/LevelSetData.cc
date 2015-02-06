
# include "../include/LevelSetData.h"

/**************************************************************************/
/*  LevelSetData.cc                                                 	  */
/*  Classe che ontiene le informazioni sul level set			 		  */
/**************************************************************************/


LevelSetData::LevelSetData ( const GetPot& dataFile,
							 const std::string& section,
							 const std::string& sectionLevelSet ):
							 M_section ( section ),
							 M_sectionLevelSet ( M_section + sectionLevelSet ),
							 M_spaceDimension ( dataFile ( ( M_section + "spaceDimension" ).data (), 1. ) ),
							 // level set
							 M_function ( dataFile ( ( M_sectionLevelSet + "levelSet" ).data (), "x" ) ),
							 M_xfunction ( dataFile ( ( M_sectionLevelSet + "xlevelSet" ).data (), "t" ) ),
							 M_yfunction ( dataFile ( ( M_sectionLevelSet + "ylevelSet" ).data (), "t" ) ),
							 M_cutFunction ( dataFile ( ( M_sectionLevelSet + "levelSetCut" ).data (), "-1" ) ),
							 M_x_map ( dataFile ( ( M_sectionLevelSet + "xMap" ).data (), "1" ) ),
							 M_y_map ( dataFile ( ( M_sectionLevelSet + "yMap" ).data (), "1" ) )
{
}



scalar_type LevelSetData::levelSetFunction ( const base_node& t )
{
    M_parser.setString ( M_function );
    M_parser.setVariable ( "x", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );

    return M_parser.evaluate ();
}// levelSetFunction


scalar_type LevelSetData::ylevelSetFunction ( const base_node& t )
{
    M_parser.setString ( M_yfunction );
    M_parser.setVariable ( "t", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );

    return M_parser.evaluate ();
}// ylevelSetFunction


scalar_type LevelSetData::levelSetCutFunction ( const base_node& t )
{
    M_parser.setString ( M_cutFunction );
    M_parser.setVariable ( "t", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );

    return M_parser.evaluate ();
}// levelSetCutFunction


scalar_type LevelSetData::y_map ( const base_node& t )
{
    M_parser.setString ( M_y_map );
    M_parser.setVariable ( "t", t [ 0 ] );

    return M_parser.evaluate ();
}// y_map


scalar_type LevelSetData::x_map ( const base_node& t )
{
    M_parser.setString ( M_x_map );
    M_parser.setVariable ( "t", t [ 0 ] );

    return M_parser.evaluate ();
}// x_map
