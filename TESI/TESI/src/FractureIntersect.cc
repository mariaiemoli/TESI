
# include "../include/FractureIntersect.h"

/**************************************************************************/
/*  FractureIntersect.cc												  */
/*  Classe che contiene tutte le intersezioni       					  */
/**************************************************************************/

/*
void FractureIntersect::constructIntesection ( const FracturePtrContainer_Type& fractures )
{
	size_type numberFractures = fractures.size();

	sizeVectorContainer_Type nodes ( numberFractures );
    sizeVectorContainer_Type listOfFractures ( numberFractures );


    //per ogni frattura cerco le intersezioni con le altre e aggiungo l'informazione che si intersecano

	for ( size_type f = 0; f < numberFractures; f++ )
	{
		// CORREGGERE!!!! in questo modo non riesco a gestire il caso in cui una frattura abbia intersezione in entrambi i nodi
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
				//fracturesInvolved[ f ]->setFractureIntersection ( nodes [ f ], *fracturesInvolved[ i ] );	// ???????????
				fracturesInvolved[ f ]->setFractureIntersection ( nodes [ f ], *fracturesInvolved[ i ] );
			}
		}
	}


	// ora che so quali sono le fratture che si intersecano e dove posso costruire le intersezioni


	// QUESTA PARTE È DA RIFARE COMPLETAMENTE!!!!!!!!!!!!!!!!!!!!!!!!

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

			// Qui faccio costruire anche il volume di intersezione
			intersection.setIntersection ( nodes[ i ], fracturesInvolved );

			std::cout << " i " << i << " nodes " << nodes [ i ] << std::endl;

			push_back ( intersection );

		}
	}

	return;

}// constructIntesection
*/

/* DOMENICA MATTINA

void FractureIntersect::constructIntesection ( const FracturePtrContainer_Type& fractures )
{
	sizeVectorContainer_Type listOfFractures ( fractures. size() );
	sizeVectorContainer_Type listOfNodes ( fractures. size() );

	for ( size_type f = 0; f < fractures.size(); f++ )
	{
		listOfFractures [ f ].clear();

		listOfNodes [ f ].clear();
		listOfNodes [ f ].resize( fractures.size () );

		// cerco tutte le possibili intersezioni tra la frattura corrente e le altre fratture
		findIntersection ( fractures [ f ], fractures, listOfNodes [ f ],  listOfFractures [ f ] );

		// ACTUNG!!! MANCA LA PARTE DI CHIAMATA ALLA FUNZIONE setFractureIntersection

	}

	//sizeVectorContainer_Type fracturesInvolved ( fractures.size() );
	std::vector<FracturePtrContainer_Type> fracturesInvolved ( fractures.size() );

	for ( size_type i = 0; i < fractures.size(); i++ )
	{
		// Prendo i puntatori alle fratture coinvolte
		//FracturePtrContainer_Type fracturesInvolved ( listOfFractures[ i ].size() +1 );

		std::cout << " fractures.size()  " << fractures.size() << std::endl;
		std::cout << " listOfFractures [ " << i << " ].size ()  " << listOfFractures [ i ].size () << std::endl;
		if ( listOfFractures [ i ].size () != 0 )	// cioè c'è un'intersezione
		{
			fracturesInvolved [ i ].clear();
			fracturesInvolved [ i ].resize( listOfFractures[ i ].size() +1 );

			fracturesInvolved [ i ][ 0 ] = fractures[ i ];

			for ( size_type f = 1; f < fracturesInvolved.size(); ++f )
			{
				 fracturesInvolved [ i ][ f ] = fractures [ listOfFractures [ i ] [ f-1 ] ];
			}

			// Costruisco la classe IntersectData per la nuova intersezione
			IntersectData intersection;

			// Qui faccio costruire anche il volume di intersezione
			intersection.setIntersection ( listOfNodes[ i ], fracturesInvolved [ i ] );

			std::cout << " listOfNodes " << listOfNodes [ i ][ 0 ] << "  " << listOfNodes [ i ][ 1 ]<< "   " << listOfNodes [ i ][ 2 ] << std::endl;

			push_back ( intersection );

			// tolgo le fratture già analizzate dagli altri vettori
			clear ( listOfFractures, fracturesInvolved [ i ] );
		}

	}


	return;
}// constructIntesection

*/


void FractureIntersect::constructIntesection ( const FracturePtrContainer_Type& fractures )
{
	sizeVectorContainer_Type listOfFractures ( fractures. size() );
	sizeVectorContainer_Type listOfNodes ( fractures. size() );

	for ( size_type f = 0; f < fractures.size(); f++ )
	{
		listOfFractures [ f ].clear();

		listOfNodes [ f ].clear();
		listOfNodes [ f ].resize( fractures.size () );

		// cerco tutte le possibili intersezioni tra la frattura corrente e le altre fratture
		findIntersection ( fractures [ f ], fractures,/* listOfNodes [ f ],*/  listOfFractures [ f ] );

		// ACTUNG!!! MANCA LA PARTE DI CHIAMATA ALLA FUNZIONE setFractureIntersection

	}

	//sizeVectorContainer_Type fracturesInvolved ( fractures.size() );
	std::vector<FracturePtrContainer_Type> fracturesInvolved ( fractures.size() );

	for ( size_type i = 0; i < fractures.size(); i++ )
	{
		// Prendo i puntatori alle fratture coinvolte
		//FracturePtrContainer_Type fracturesInvolved ( listOfFractures[ i ].size() +1 );

		if ( listOfFractures [ i ].size () != 0 )	// cioè c'è un'intersezione
		{
			fracturesInvolved [ i ].clear();
			fracturesInvolved [ i ].resize( listOfFractures[ i ].size() +1 );

			fracturesInvolved [ i ][ 0 ] = fractures[ i ];

			for ( size_type f = 1; f < fracturesInvolved.size(); ++f )
			{
				 fracturesInvolved [ i ][ f ] = fractures [ listOfFractures [ i ] [ f-1 ] ];
			}
		}
	}

	for ( size_type i = 0; i < fractures.size(); i++ )
	{

		// Costruisco la classe IntersectData per la nuova intersezione
		IntersectData intersection;

		// Qui faccio costruire anche il volume di intersezione
		intersection.setIntersection (/* listOfNodes[ i ], */fracturesInvolved [ i ] );

		push_back ( intersection );

		// tolgo le fratture già analizzate dagli altri vettori
		clear ( listOfFractures, fracturesInvolved [ i ] );
	}



	return;
}// constructIntesection


void FractureIntersect::findIntersection ( const FractureHandlerPtr_Type& f, const FracturePtrContainer_Type& fractures,
										   /*sizeVector_Type& listOfNodes, */ sizeVector_Type& listOfFractures )
{
	size_type numberFracture = fractures.size();
	size_type ID = f->getID();
	size_type nbDof = f->getData().getSpatialDiscretization ()-1;
//	nodes.clear();

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

			/*
			std::cout << " otherFracture " << otherFracture << std::endl;
			std::cout << " node_of[ 0 ] " << node_of[ 0 ][0]  << "   " << node_of[ 0 ][ 1 ] << std::endl;
			std::cout << " node_of[ 1 ] " << node_of[ 1 ][0]  << "   " << node_of[ 1 ][ 1 ] << std::endl;
			 */

			if( node_of[ 0 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				std::cout << " caso 1  a " << std::endl;
				//listOfNodes [ otherFracture ] = 0;	// metto in coda il nodo di f in cui c'è l'intersezione
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				std::cout << " caso 2  a " << std::endl;
				//listOfNodes [ otherFracture ] = 0;
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 0 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 0 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				std::cout << " caso 1  b " << std::endl;
				//listOfNodes [ otherFracture ] = nbDof-1;
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 1 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 1 ][ 1 ] )
			{
				std::cout << " caso 2  b " << std::endl;
				//listOfNodes [ otherFracture ] = nbDof-1;
				listOfFractures.push_back( otherFracture );
			}

		}

	}

	return;

}// findIntersection

/*
void FractureIntersect::findIntersection ( const FractureHandlerPtr_Type& f, const FracturePtrContainer_Type& fractures,
										   sizeVector_Type& listOfNodes, sizeVector_Type& listOfFractures )
{
	size_type numberFracture = fractures.size();
	size_type ID = f->getID();
	size_type nbDof = f->getData().getSpatialDiscretization ()-1;
//	nodes.clear();

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
				listOfNodes [ otherFracture ] = 0;	// metto in coda il nodo di f in cui c'è l'intersezione
				listOfFractures.push_back( otherFracture );
			}
			else if( node_of[ 1 ][ 0 ] == node_f[ 0 ][ 0 ] && node_of[ 1 ][ 1 ] == node_f[ 0 ][ 1 ] )
			{
				listOfNodes [ otherFracture ] = 0;
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

	return;

}// findIntersection
*/

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
			if( !M_intersections[ i ].isEqual( intersection ))
			{
				/*M_intersections[ i ].updateNodes( intersection );
			}
			else
			{*/
				M_intersections.push_back( intersection );
			}
		}
	}

	return;
}// push_back


void FractureIntersect::clear ( sizeVectorContainer_Type& listOfFractures, const FracturePtrContainer_Type& fracturesInvolved )
{
	sizeVector_Type fractures_id;

	for ( size_type i = 0; i < fracturesInvolved.size(); i++ )
	{
		fractures_id.push_back( fracturesInvolved [ i ]->getID () );
	}

	for( size_type f = 0; f < listOfFractures.size(); f++ )
	{
		for ( size_type i = 0; i < fractures_id.size(); i++ )
		{
			for ( size_type j = 0; j < listOfFractures [ f ].size(); j++ )
			{
				if ( listOfFractures [ f ][ j ] == fractures_id [ i ] )
				{
					listOfFractures [ f ].erase ( listOfFractures [ f ].begin() + j );
				}
			}
		}
	}

	return;
}// clear


IntersectDataContainer_Type FractureIntersect::getIntersections () const
{
    return M_intersections;

} // getIntersections


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


size_type FractureIntersect::getNumberIntersection() const
{
	return M_intersections.size();
}// getNumberIntersection



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
