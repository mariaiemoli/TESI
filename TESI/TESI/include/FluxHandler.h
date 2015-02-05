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

#ifndef _FLUXHANDLER_H_
#define _FLUXHANDLER_H_ 1

/**************************************************************************/
/*  FluxHandler.h                                                   	  */
/*  Classe che rappresenta una funzione flusso					 		  */
/**************************************************************************/


#include "StringUtility.h"
#include "Core.h"
#include "Parser.h"


class FluxHandler
{
public:

	FluxHandler ( const GetPot& dataFile, const std::string& sectionFlux );


	/**
	 * Funzione che aggiorna il valore della condizone al bordo nel nodo di intersezione
	 */
	void update_UI ( const size_type& i, const scalar_type& u );


	inline std::string getFlux() const
	{
		return M_flux;
	}

	inline std::string getFlux1() const
	{
		return M_flux1;
	}

	inline std::string getMonotone() const
	{
		return M_Monotone;
	}

	inline scalar_type getH() const
	{
		return M_H;
	}

	inline scalar_type getUs() const
	{
		return M_Us;
	}

	inline scalar_type getFirst() const
	{
		return M_first;
	}

	inline scalar_type getSecond() const
	{
		return M_second;
	}

	inline scalar_type geta() const
	{
		return M_a;
	}

	inline scalar_type getb() const
	{
		return M_b;
	}


private:

	std::string M_sectionFlux;

	std::string M_flux;
	std::string M_flux1;

	std::string M_Monotone;

	size_type M_H;
	scalar_type M_Us;
	scalar_type M_first;
	scalar_type M_second;

	scalar_type M_a;
	scalar_type M_b;

	mutable LifeV::Parser M_parser;
};





typedef FluxHandler FluxHandler_Type;												/*!< Classe FluxHandler */
typedef boost::shared_ptr < FluxHandler_Type > FluxHandlerPtr_Type;					/*!< Puntatore alla classe FluxHandler */
typedef std::vector<FluxHandlerPtr_Type> FluxPtrContainer_Type;						/*!< vettore di puntatori alla classe FluxHandler */
typedef std::vector<FluxHandlerPtr_Type> FluxContainer_Type;						/*!< vettore di classi FluxHandler */
typedef boost::shared_ptr<FluxPtrContainer_Type> FluxPtrContainerPtr_Type;			/*!< puntatore a un vettore di puntatori alla classe FluxHandler */



#endif /* _FLUXHANDLER_H_ */
