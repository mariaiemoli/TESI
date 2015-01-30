/* Progetto di tesi
 * Politecnico di Milano
 *
 * \author Iemoli Maria
 * \date 2015
 *
 * Problema di saturazione per un fluido bifase in un network di fratture
 *
 */

#ifndef _FRACTUREINTERSECT_
#define _FRACTUREINTERSECT_ 1

#include "Core.h"
#include "IntersectData.h"
#include "FractureHandler.h"
#include "UsefulFunctions.h"
#include <assert.h>
#include <map>

/**************************************************************************/
/*  FractureIntersect.h													  */
/*  Classe che contiene tutte le intersezioni       					  */
/**************************************************************************/

class FractureIntersect
{
public:

	FractureIntersect ()
	{};

	/**
	 * Con questa funzione costruisco la classe che contiene tutte le intersezioni valutando, per ogni frattura singolarmente,
	 * se esiste un punto in comune con le altre.
	 * Poichè le fratture sono definite come segmenti non tagliati, le intersezioni tra fratture possono trovarsi solo a inizio
	 * o fine frattura
	 */
	void constructIntesection ( const FracturePtrContainer_Type& fractures );

	/**
	 * Con questa funzione verifico se una data frattura ha intersezione con tutte le altre.
	 * Quello che ottengo è un vettore di nodi che indica dove la frattura ha intersezione (per ora sarà un nodo solo
	 * e, in particolare, sarà il primo o l'ultimo) e la lista degli indice delle fratture con cui ha intersezione
	 */
	void findIntersection ( const FractureHandlerPtr_Type& f, const FracturePtrContainer_Type& fractures,
							VectorIntersection_Type& listOfFractures );


	/**
	 * Dopo aver costruito una nuova intersezione, prima di inserirla nel vettore delle intersezioni, verifico che non vi sia già.
	 * In tal caso non inserisco le fratture coinvolte, ma solo gli indici dei nodi corrispondenti non ancora inseriti
	 */
	void Push_Back( const IntersectData& intersection );

	/**
	 * Funzione che elimina dalla lista delle fratture con intersezioni quelle già analizzate
	 */
	void clear ( VectorIntersectionContainer_Type& listOfFractures, const size_type& i, const size_type& j );


	bool isIn ( const size_type& i, const sizeVector_Type& fractures) const;


    /**
     * Funzione che restituisce tutte le intersezioni.
     * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni
     */
    IntersectDataContainer_Type getIntersections () const
    {
        return M_intersections;
    }


    /**
     * Funzione che restituisce tutte le intersezioni di tipo " Cross ".
     * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni del tipo " Cross "
     */
    IntersectDataContainer_Type getCrossIntersections () const;


    /**
     * Funzione che restituisce tutte le intersezioni di tipo " Bifurcation ".
     * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni del tipo " Bifurcation "
     */
    IntersectDataContainer_Type getBifurcationIntersections () const;

	size_type getNumberIntersection () const
	{
		return M_intersections.size();
	}

	size_type getNumberCross () const;

	size_type getNumberBifurcation () const;


private:

	IntersectDataContainer_Type M_intersections;

};

typedef FractureIntersect FractureIntersect_Type;									/*!< classe FractureIntersect */
typedef boost::shared_ptr < FractureIntersect > FractureIntersectPtr_Type;			/*!< puntatore alla classe FractureIntersect */

#endif /* _FRACTUREINTERSECT_H_ */
