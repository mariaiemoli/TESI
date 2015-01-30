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



#ifndef _INTERSECTDATA_
#define _INTERSECTDATA_ 1

#include "FractureHandler.h"
#include "Geometry.h"

/**************************************************************************/
/*  IntersectData.h														  */
/*  Classe che contiene tutte le informazioni su un'intersezione		  */
/**************************************************************************/

class IntersectData
{
public:

	IntersectData ()
	{};


    IntersectData ( const IntersectData& in )
    {
    	copy ( in );
    } // costruttore di copia


    /**
     * Funzione che definisce l'intersezione: definisce le fratture coinvolte e assegna al vettore dei nodi la dimensione
     * corretta, ossia pari al numero di fratture. Questo perchè per ogni frattura coinvolta nell'intersezione voglio salvare il
     * corrispondente grado di libertà dove avviene l'incontro.
     */
    void setIntersection ( const FracturePtrContainer_Type& fractures );

    /**
     * Funzione che verifica se due intersezioni sono uguali a meno dei nodi, ossia se le fratture coinvolte sono le stesse
     */
    bool isEqual ( const IntersectData& intersection );

    /**
     * Funzione che aggiunge uno degli indici dei nodi di intersezione mancanti
     */
    void updateNodes ();

    /**
     * Funzione che aggiorna le bc per le fratture che si intersecano
     */
    void update_UI ( const size_type& i );


    /**
     * Funzione che aggiorna il valore della saturazione all'intersezione
     */
    void update_Ui ( const scalar_type& landa);

    const FracturePtrContainer_Type& getFractures () const
    {
        return M_fractures;
    } // getFractures

    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
        return M_fractures [ f ];
    } // getFracture


    scalar_type measure () const
    {
        return M_measure;
    } // measure

    size_type getNumFractures () const
    {
        return M_fractures.size();
    } // getNumFractures

    sizeVector_Type getIntersectionPoint () const
    {
    	return M_intersectionPoint;
    }


    scalar_type getU0 () const
    {
    	return M_u0;
    }


    const Intersection_Type& getIntersect () const
    {
    	return M_intersection;
    }

	//Operatori
    IntersectData& operator=(const IntersectData& t);


private:

    void copy ( const IntersectData& in );

    sizeVector_Type M_intersectionPoint;

    // Insieme delle fratture che si intersecano
    FracturePtrContainer_Type M_fractures;

    scalar_type M_u0;

    scalar_type M_measure;

	Intersection_Type M_intersection;

};

typedef IntersectData IntersectData_Type;											/*!< Classe IntersectData */
typedef std::vector < IntersectData_Type > IntersectDataContainer_Type;				/*!< Vettore di classi IntersectData */
typedef boost::shared_ptr < IntersectData_Type > IntersectDataPtr_Type;				/*!< Puntatore alla classe IntersectData */
typedef std::vector < IntersectDataPtr_Type > IntersectDataPtrContainer_Type;		/*!< Vettore di puntatori alla classe IntersectData */


#endif /* _INTERSECTDATA_H_ */
