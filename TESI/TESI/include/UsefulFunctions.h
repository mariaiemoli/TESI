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

/**
 * Funzioni per esportare soluzioni e mesh
 */
void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const scalarVector_Type& solution );


void exportMesh ( const std::string& fileName, const getfem::mesh& mesh );

#endif /* _USEFULFUNCTION_H_ */