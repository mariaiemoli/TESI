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

#ifndef _BCHANDLER_H_
#define _BCHANDLER_H_ 1

/**************************************************************************/
/*  BCHandler.h          		                                          */
/*  Classe in cui vengono imposte le condizioni al contorno				  */
/**************************************************************************/

#include "Core.h"
#include "FractureData.h"
#include "BC.h"



class BCHandler
{
public:

    // Costruttore
    BCHandler ( const BCPtrContainer_Type& fractureBC );


    inline const BCPtr_Type& getFractureBC ( const size_type& id ) const
    {
        return M_fractureBC [ id ];
    }


private:

	 BCPtrContainer_Type M_fractureBC;
};

typedef BCHandler BCHandler_Type;											/*!< Classe BCHandler */
typedef boost::shared_ptr<BCHandler_Type> BCHandlerPtr_Type;				/*!< Puntatore alla classe BCHandler */

#endif /* BCHANDLER_H_ */
