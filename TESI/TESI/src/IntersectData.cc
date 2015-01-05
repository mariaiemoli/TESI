
# include "../include/IntersectData.h"

/**************************************************************************/
/*  IntersectData.cc													  */
/*  Classe che contiene tutte le informazioni su un'intersezione		  */
/**************************************************************************/

void IntersectData::copy ( const IntersectData& in )
{
    if ( this != &in )
    {
        M_fractures.resize ( in.M_fractures.size() );
        for ( size_type i = 0; i < M_fractures.size(); ++i )
        {
            M_fractures[i] = in.M_fractures[i];
        }

        M_intersectionPoint = in.M_intersectionPoint;
    }

    return;
} // copy


void IntersectData::setIntersection ( const sizeVector_Type& nodes,
					   	   	   	   	  const FracturePtrContainer_Type& fractures )
{
	M_fractures = fractures;
	M_intersectionPoint.clear();
	M_intersectionPoint.resize( fractures.size() );

	for ( size_type i = 0; i < nodes.size(); i++ )
	{
		M_intersectionPoint[ i ] = nodes[ i ];
	}

	return;
}// setIntersection


bool IntersectData::isEqual ( const IntersectData& intersection )
{
	size_type count = 0;

	if ( M_fractures.size() != intersection.getFractures ().size() )
	{
		return false;
	}
	else
	{
		for ( size_type i = 0; i < M_fractures.size(); i++ )
		{
			for ( size_type j = 0; j < intersection.getFractures ().size(); j++ )
			{
				if ( M_fractures[ i ]->getID() == intersection.getFracture ( j )->getID () )
				{
					count ++;
				}
			}
		}

		return ( count == M_fractures.size() );
	}
}// isEqual


void IntersectData::updateNodes ( const IntersectData& intersection )
{
	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		for ( size_type j = 0; j < intersection.getFractures ().size(); j++ )
		{
			if ( M_fractures[ i ]->getID() == intersection.getFracture ( j )->getID () )
			{
				if ( intersection.getIntersectionPoint ().size() != 0 )
				{
					M_intersectionPoint [ i ] = intersection.getIntersectionPoint () [ j ];
				}
			}
		}
	}

	return;
}// updateNodes
