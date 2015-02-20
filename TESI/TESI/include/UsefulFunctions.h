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

#ifndef _USEFULFUNCTION_H_
#define _USEFULFUNCTION_H_ 1

/**************************************************************************/
/*  UsefulFunction.h                                                 	  */
/*  Libreria in cui definisco funzioni ausiliarie						  */
/**************************************************************************/

#include "Core.h"
#include "Parser.h"


/**
 * Funzioni per esportare soluzioni e mesh
 */
void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution );


void exportMesh ( const std::string& fileName, const getfem::mesh& mesh );

scalar_type min ( const std::string& flusso, const scalar_type ul, const scalar_type ur );

scalar_type max ( const std::string& flusso, const scalar_type ul, const scalar_type ur );

/**
 * Funzione che calcola la distanza tra due punti
 */
scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 );


class VectorIntersection
{
public:

	VectorIntersection ();


	sizeVector_Type const & operator[]( size_type i ) const
	{
		if ( i == 0 )
		{
			return Vector_0;
		}
		else
		{
			return Vector_1;
		}
	}


	sizeVector_Type & operator[]( size_type i )
	{
		if ( i == 0 )
		{
			return Vector_0;
		}
		else
		{
			return Vector_1;
		}
	}



private:

	sizeVector_Type Vector_0;

	sizeVector_Type Vector_1;
};

typedef VectorIntersection VectorIntersection_Type;									/*!< classe VectorIntersection */
typedef std::vector < VectorIntersection_Type > VectorIntersectionContainer_Type;				/*!< Vettore di classi IntersectData */
typedef boost::shared_ptr < VectorIntersection > VectorIntersectionPtr_Type;			/*!< puntatore alla classe VectorIntersection */


void orderId( size_type& id_i, size_type& id_j, size_type& id_k );

#endif /* _USEFULFUNCTION_H_ */
