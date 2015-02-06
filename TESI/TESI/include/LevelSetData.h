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

#ifndef _LEVELSETDATA_H_
#define _LEVELSETDATA_H_ 1


#include "Core.h"
#include "Parser.h"


/**************************************************************************/
/*  LevelSetData.h                                                  	  */
/*  Classe che ontiene le informazioni sul level set			 		  */
/**************************************************************************/


class LevelSetData
{
public:

	// Costruttore
    LevelSetData ( const GetPot& dataFile,
                   const std::string& section, // = "fractureData/",
                   const std::string& sectionLevelSet = "levelSet/" );

    /**
     *
     * Questa funzione definisce il level set che rappresenta la frattura, level set valutato in (x,y).
     * \param base_node& x: nodo in coordinate ( x, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     *
     */
    scalar_type levelSetFunction ( const base_node& x );


    /**
     *
     * Questa funzione definisce il level set che rappresenta la frattura, level set valutato in (t,y).
     * \param base_node& x: nodo in coordinate ( t, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     */
    scalar_type ylevelSetFunction ( const base_node& x );


    /** scalar_type levelSetCutFunction ( const base_node& x, int num = 0 )
     *
     * Questa funzione definisce il  level set che rappresenta la frattura, se per caso voglio "tagliare la frattura".
     * \param base_node& x: nodo in coordinate ( t, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     */
    scalar_type levelSetCutFunction ( const base_node& x );


    /**
     * Questa funzione rappresenta la mappa dalla frattura piatta, y(t).
     */
    scalar_type y_map ( const base_node& t );


    /**
	 * Questa funzione rappresenta la mappa dalla frattura piatta, x(t).
	 */
    scalar_type x_map ( const base_node& t );


    inline bgeot::dim_type getSpaceDimension () const
    {
        return M_spaceDimension;
    }

    std::string getLevelSetFunctionString () const
    {
        return M_yfunction;
    }



private:

    // Attributes
    std::string M_section;
    std::string M_sectionLevelSet;

    bgeot::dim_type M_spaceDimension;

    std::string M_function;
    std::string M_xfunction;
    std::string M_yfunction;
    std::string M_cutFunction;

    std::string M_x_map;
    std::string M_y_map;

    mutable LifeV::Parser M_parser;

};


typedef LevelSetData LevelSetData_Type;								/*!< Classe LevelSetData */
typedef boost::shared_ptr<LevelSetData_Type> LevelSetDataPtr_Type;	/*!< Puntatore alla classe LevelSetData */


#endif /* _LEVELSETDATA_H_ */
