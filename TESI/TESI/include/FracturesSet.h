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

#ifndef _FRACTURESSET_H_
#define _FRACTURESSET_H_ 1

/**************************************************************************/
/*  FracturesSet.h                                                   	  */
/*  Classe che rappresenta l'insieme di tutte le fratture		 		  */
/**************************************************************************/

#include "Core.h"
#include "FractureHandler.h"
#include "FractureIntersect.h"

class FracturesSet
{
public:

	// Costruttore nullo
	FracturesSet ();

	void init ( const GetPot& dataFile, const std::string& section, const size_type& numFractures, const ExporterPtr_Type& exporter );

    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
            return M_fractures[f];
    }

	inline size_type getNumberFractures() const
	{
		return M_fractures.size();
	}

	const FractureIntersectPtr_Type& getIntersections() const
	{
		return M_intersections;
	}


private:

	// vettore di puntatori a tutte le fratture
	FracturePtrContainer_Type M_fractures;

	// intersezioni
	FractureIntersectPtr_Type M_intersections;


};

typedef FracturesSet FracturesSet_Type;										/*!< Classe FracturesSet */
typedef boost::shared_ptr < FracturesSet_Type > FracturesSetPtr_Type;		/*!< Puntatore alla classe FracturesSet */


#endif /* _FRACTURESSET_H_ */
