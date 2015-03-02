
# include "../include/FractureIntersect.h"

/**************************************************************************/
/*  FractureIntersect.cc												  */
/*  Classe che contiene tutte le intersezioni       					  */
/**************************************************************************/


void FractureIntersect::constructIntesection ( const FracturePtrContainer_Type& fractures )
{
	VectorIntersectionContainer_Type listOfFractures ( fractures. size() );

	for ( size_type f = 0; f < fractures.size(); f++ )
	{
		// cerco tutte le possibili intersezioni tra la frattura corrente e le altre fratture
		findIntersection ( fractures [ f ], fractures, listOfFractures [ f ] );
	}

	std::vector<FracturePtrContainer_Type> fracturesInvolved ( fractures.size() );


	for ( size_type i = 0; i < fractures.size(); i++ )
	{
		// Prendo i puntatori alle fratture coinvolte

		for ( size_type j = 0; j < 2; j++ )
		{
			if ( listOfFractures [ i ] [ j ].size () != 0 )	// cioè se la frattura i-esima ha un'intersezione nel nodo j
			{
				fracturesInvolved [ i ].clear();
				fracturesInvolved [ i ].resize( listOfFractures[ i ] [ j ].size() +1 );

				fracturesInvolved [ i ][ 0 ] = fractures[ i ];

				for ( size_type f = 1; f < fracturesInvolved.size(); ++f )
				{
					 fracturesInvolved [ i ][ f ] = fractures [ listOfFractures [ i ] [ j ] [ f-1 ] ];
				}

				// Costruisco la classe IntersectData per la nuova intersezione
				IntersectData intersection;

				// Qui faccio costruire anche il volume di intersezione
				intersection.setIntersection ( fracturesInvolved [ i ] );

				Push_Back ( intersection );

				// tolgo le fratture già analizzate dagli altri vettori
				clear ( listOfFractures,  i, j );

			}
		}
	}

	// ************************************ //


	listOfFractures.clear();
	listOfFractures.resize( fractures. size() );

	for ( size_type f = 0; f < fractures.size(); f++ )
	{
		// cerco tutte le possibili intersezioni tra la frattura corrente e le altre fratture
		findIntersection ( fractures [ f ], fractures, listOfFractures [ f ] );
	}

	for ( size_type i = 0; i < fractures.size(); i++ )
	{
		// Prendo i puntatori alle fratture coinvolte
		for ( size_type j = 0; j < 2; j++ )
		{
			if ( listOfFractures [ i ] [ j ].size () != 0 )	// cioè se la frattura i-esima ha un'intersezione nel nodo j
			{
				fracturesInvolved [ i ].clear();
				fracturesInvolved [ i ].resize( listOfFractures[ i ] [ j ].size() +1 );

				fracturesInvolved [ i ][ 0 ] = fractures[ i ];

				for ( size_type f = 1; f < fracturesInvolved.size(); ++f )
				{
					 fracturesInvolved [ i ][ f ] = fractures [ listOfFractures [ i ] [ j ] [ f-1 ] ];
				}

				for ( size_type k = 1; k < fracturesInvolved[ i ].size(); k++ )
				{
					fractures [ i ]->setMeshLevelSetFracture ( *fracturesInvolved [ i ][ k ] );
				}

			}
		}
	}

	return;
}// constructIntesection



void FractureIntersect::findIntersection ( const FractureHandlerPtr_Type& f,
										   const FracturePtrContainer_Type& fractures,
										   VectorIntersection_Type& listOfFractures )
{
	size_type numberFracture = fractures.size();
	size_type ID = f->getID();

	scalarVectorContainer_Type node_f( 2 );
	scalarVectorContainer_Type node_of( 2 );
	base_node t0(1);
	base_node t1(1);
	t0[ 0 ] = 0;
	t1[ 0 ] = 1;

	node_f[ 0 ].push_back( f->getData ().getA ());
	node_f[ 1 ].push_back( f->getData ().getB ());
	node_f[ 0 ].push_back( f->getLevelSet ()->getData ( )->y_map( t0 ));
	node_f[ 1 ].push_back( f->getLevelSet ()->getData ( )->y_map( t1 ));

	for ( size_type otherFracture = 0; otherFracture < numberFracture; otherFracture++ )
	{
		node_of.clear();
		node_of.resize( 2 );


		if( otherFracture != ID )
		{
			node_of[ 0 ].push_back( fractures[ otherFracture ]->getData ().getA ());
			node_of[ 1 ].push_back( fractures[ otherFracture ]->getData ().getB ());
			node_of[ 0 ].push_back( fractures[ otherFracture ]->getLevelSet ()->getData ( )->y_map( t0 ));
			node_of[ 1 ].push_back( fractures[ otherFracture ]->getLevelSet ()->getData ( )->y_map( t1 ));


			if( node_of[ 0 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				listOfFractures [ 0 ].push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				listOfFractures [ 0 ].push_back( otherFracture );
			}
			else if( node_of[ 0 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				listOfFractures [ 1 ].push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				listOfFractures [ 1 ].push_back( otherFracture );
			}

		}

	}

	return;

}// findIntersection



void FractureIntersect::Push_Back( const IntersectData& intersection )
{
	if ( M_intersections.size() == 0 )
	{
		IntersectData tmp ( intersection );

		M_intersections.push_back( tmp );

	}
	else
	{
		for ( size_type i = 0; i < M_intersections.size(); i++ )
		{
			if( !M_intersections[ i ].isEqual( intersection ))
			{
				IntersectData tmp = intersection;

				M_intersections.push_back( tmp );
			}
		}
	}

	return;
}// push_back


void FractureIntersect::clear ( VectorIntersectionContainer_Type& listOfFractures, const size_type& i, const size_type& j )
{
	sizeVector_Type fractures = listOfFractures [ i ][ j ];

	for ( size_type k = 0; k < listOfFractures.size(); k++ )
	{
		if ( k!= i )
		{
			if ( isIn ( k, fractures ))
			{
				if ( isIn( i, listOfFractures [ k ][ 0 ]))
				{
					listOfFractures [ k ][ 0 ].clear();
				}
				else
				{
					listOfFractures [ k ][ 1 ].clear();
				}
			}
		}
	}

	return;
}// clear



bool FractureIntersect::isIn ( const size_type& i, const sizeVector_Type& fractures) const
{
	for ( size_type j = 0; j < fractures.size(); j++ )
	{
		if ( fractures [ j ] == i )
		{
			return true;
		}
	}

	return false;
}// isIn



IntersectDataContainer_Type FractureIntersect::getCrossIntersections () const
{
	IntersectDataContainer_Type tmp;
	size_type Num_Cross = getNumberCross();

	if( Num_Cross !=0 )
	{
		for ( size_type i = 0; i< M_intersections.size(); i++ )
		{
			if ( M_intersections[ i ].getNumFractures () == 4 )
			{
				tmp.push_back( M_intersections[ i ] );
			}
		}
	}

    return tmp;

} // getCrossIntersections


IntersectDataContainer_Type FractureIntersect::getBifurcationIntersections () const
{
	IntersectDataContainer_Type tmp;

	size_type Num_Bifu = getNumberBifurcation();

	if( Num_Bifu !=0 )
	{
		for ( size_type i = 0; i< M_intersections.size(); i++ )
		{
			if ( M_intersections[ i ].getNumFractures () == 3 )
			{
				tmp.push_back( M_intersections[ i ] );
			}
		}
	}

    return tmp;

} // getBifurcationIntersections


size_type FractureIntersect::getNumberCross() const
{
	size_type count = 0;

	for ( size_type i = 0; i < M_intersections.size (); i++ )
	{
		if ( M_intersections[ i ].getFractures().size() == 4 )
		{
			count ++;
		}
	}

	return count;
}// getNumberCross

size_type FractureIntersect::getNumberBifurcation() const
{
	size_type count = 0;

	for ( size_type i = 0; i < M_intersections.size (); i++ )
	{
		if ( M_intersections[ i ].getFractures().size() == 3 )
		{
			count ++;
		}
	}

	return count;
}// getNumberBifurcation
