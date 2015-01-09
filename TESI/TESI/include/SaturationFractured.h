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

#ifndef _SATURATIONFRACTURED_H_
#define _SATURATIONFRACTURED_H_ 1

/**************************************************************************/
/*  SaturationFractured.h                                                 */
/*  Classe in cui viene risolto il problema di saturazione				  */
/**************************************************************************/

#include "Core.h"
#include "FracturesSet.h"
#include "StringUtility.h"
#include "Parser.h"
#include "Exporter.h"
#include "BCHandler.h"
#include "UsefulFunctions.h"
#include <algorithm>
#include <math.h>


class SaturationFractured
{
public:
	SaturationFractured( const GetPot& dataFile,
						 FracturesSetPtr_Type& fractures,
						 const BCHandlerPtr_Type& bcHandler,
						 const ExporterPtr_Type& exporter,
					     const std::string time = "dataTime/" );

	/**
	 * Funzione che inizializza le matrici per costruire il sistema algebrico
	 */
	void init();


	/**
	 * Funzione che assembla il sistema algebrico
	 */
	void assembly( const scalar_type& landa, const scalarVectorContainer_Type& u0, const scalarVectorContainer_Type& flux );


	void solve();


	/**
	 * Funzione che risolve il problema nel caso un cui nella frattura in questione non vi siano discontinuità nel materiale
	 */
	void solve_continuity ( const size_type f, const scalarVector_Type& u0, scalarVector_Type& Flux, const size_type k = 0 );


	/**
	 * Funzione che risolve il problema nel caso un cui nella frattura in questione vi siano discontinuità nel materiale
	 */
	void solve_discontinuity ( const size_type f,  const scalarVector_Type& u0, scalarVector_Type& Flux );


	FracturesSetPtr_Type& get_Fracture()
	{
		return M_fractures;
	}

	inline scalar_type getFinalTime () const
	{
		return M_t;
	}

private:

	std::string M_section;

	FracturesSetPtr_Type M_fractures;

    // BC Handler
    BCHandlerPtr_Type M_bcHandler;

    // Exporter
    ExporterPtr_Type M_exporter;

	scalar_type M_t;

    // Global matrix, saturation
    sparseMatrixPtr_Type M_globalMatrix;
    // Termine noto di destra del sistema
    scalarVectorPtr_Type M_globalRightHandSide;
    // Soluzione del sistema
    scalarVectorPtr_Type M_saturation;

    // Soluzione del sistema nelle fratture
    scalarVectorPtrContainer_Type M_fractureSaturation;

    mutable LifeV::Parser M_parser;

};

typedef SaturationFractured SaturationFractured_Type;											/*!< Classe SaturationFractured */
typedef boost::shared_ptr<SaturationFractured_Type> SaturationFracturedPtr_Type;				/*!< Puntatore alla classe SaturationFractured */


#endif /* _SATURATIONFRACTURED_H_ */
