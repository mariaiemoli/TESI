
# include "../include/FractureIntersect.h"

/**************************************************************************/
/*  FractureIntersect.cc												  */
/*  Classe che contiene tutte le intersezioni       					  */
/**************************************************************************/

void FractureIntersect::constructIntesection ( const FracturePtrContainer_Type& fractures )
{
	size_type numberFractures = fractures.size();

	sizeVectorContainer_Type nodes ( numberFractures );
    sizeVectorContainer_Type listOfFractures ( numberFractures );

    /*
     * per ogni frattura cerco le intersezioni con le altre e aggiungo l'informazione che si intersecano
     */
	for ( size_type f = 0; f < numberFractures; f++ )
	{
		findIntersection ( fractures [ f ], fractures, nodes [ f ],  listOfFractures [ f ] );

		FracturePtrContainer_Type fracturesInvolved ( listOfFractures[ f ].size() +1 );

		fracturesInvolved[ 0 ] = fractures [ f ];

		for ( size_type j = 1; j < fracturesInvolved.size(); ++j )
		{
			 fracturesInvolved [ j ] = fractures [ listOfFractures [ f ] [ j-1 ] ];
		}

		if ( nodes[ f ].size() != 0 )
		{
			for ( size_type i = 1; i < fracturesInvolved.size(); i++ )
			{
				fracturesInvolved[ f ]->setFractureIntersection ( nodes [ f ], *fracturesInvolved[ i ] );
			}
		}
	}

	/*
	 * ora che so quali sono le fratture che si intersecano e dove posso costruire le intersezioni
	 */
	for ( size_type i = 0; i < nodes.size(); i++ )
	{
		if ( nodes [ i ].size() != 0 )
		{
			// Prendo i puntatori alle fratture coinvolte
			FracturePtrContainer_Type fracturesInvolved ( listOfFractures[ i ].size() +1 );

			fracturesInvolved[ 0 ] = fractures[ i ];

			for ( size_type f = 1; f < fracturesInvolved.size(); ++f )
			{
				 fracturesInvolved [ f ] = fractures [ listOfFractures [ i ] [ f-1 ] ];
			}

			// Costruisco la classe IntersectData per la nuova intersezione
			IntersectData intersection;

			intersection.setIntersection ( nodes[ i ], fracturesInvolved );

			push_back ( intersection );

		}
	}


}// constructIntesection


void FractureIntersect::findIntersection ( const FractureHandlerPtr_Type& f, const FracturePtrContainer_Type& fractures,
										   sizeVector_Type& nodes, sizeVector_Type& listOfFractures )
{
	size_type numberFracture = fractures.size();
	size_type ID = f->getID();
	size_type nbDof = f->getData().getSpatialDiscretization ();
	nodes.clear();

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

		if( otherFracture != ID )
		{
			node_of[ 0 ].push_back( fractures[ otherFracture ]->getData ().getA ());
			node_of[ 1 ].push_back( fractures[ otherFracture ]->getData ().getB ());
			node_of[ 0 ].push_back( fractures[ otherFracture ]->getLevelSet ()->getData ( )->y_map( t0 ));
			node_of[ 1 ].push_back( fractures[ otherFracture ]->getLevelSet ()->getData ( )->y_map( t1 ));

			if( node_of[ 0 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				nodes.push_back( 0 );
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				nodes.push_back( 0 );
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 0 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				nodes.push_back( nbDof-1 );
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				nodes.push_back( nbDof-1 );
				listOfFractures.push_back( otherFracture );
			}



		}

	}

}// findIntersection


void FractureIntersect::push_back( const IntersectData& intersection )
{
	if ( M_intersections.size() == 0 )
	{
		M_intersections.push_back( intersection );
	}
	else
	{
		for ( size_type i = 0; i < M_intersections.size(); i++ )
		{
			if( M_intersections[ i ].isEqual( intersection ))
			{
				M_intersections[ i ].updateNodes( intersection );
			}
			else
			{
				M_intersections.push_back( intersection );
			}
		}
	}

	return;
}// push_back



IntersectDataContainer_Type FractureIntersect::getCrossIntersections () const
{
	IntersectDataContainer_Type tmp;
	size_type Num_Cross = getNumberCross();

	if( Num_Cross !=0 )
	{
		for ( size_type i = 0; i< M_intersections.size(); i++ )
		{
			if ( M_intersections[ i ].getNumFractures () == 2 )
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
		if ( M_intersections[ i ].getFractures().size() == 2 )
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
