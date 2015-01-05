
# include "../include/BCHandler.h"

/**************************************************************************/
/*  BCHandler.cc          		                                          */
/*  Classe in cui vengono imposte le condizioni al contorno				  */
/**************************************************************************/

BCHandler::BCHandler ( const BCPtrContainer_Type& fractureBC ) :
                       M_fractureBC(fractureBC)
{
}
