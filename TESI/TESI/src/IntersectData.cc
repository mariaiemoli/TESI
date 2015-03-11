
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

        M_measure = in.M_measure;

        M_u0 = in.M_u0;

        M_intersection = in.M_intersection;
    }

    return;
} // copy


void IntersectData::setIntersection ( const FracturePtrContainer_Type& fractures )
{
	M_fractures = fractures;
	M_intersectionPoint.clear();
	M_intersectionPoint.resize( fractures.size() );

	updateNodes();

	// devo aggiornare anche le bc in questo caso

	FracturePtrContainer_Type fractures_copy = fractures;

	M_intersection.setIntersection ( fractures_copy );

	M_u0 = 0.;

	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		M_fractures [ i ] -> getData().getFluxHandler( 0 )->update_UI ( M_intersectionPoint [ i ], M_u0 );
	}

	M_measure = M_intersection.measure();

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


void IntersectData::updateNodes ()
{
	if ( M_fractures.size() == 3 )
	{
		M_intersection.Resize(3);

		size_type nbDof0 = M_fractures [ 0 ]->getData().getSpatialDiscretization ()-1;
		size_type nbDof1 = M_fractures [ 1 ]->getData().getSpatialDiscretization ()-1;
		size_type nbDof2 = M_fractures [ 2 ]->getData().getSpatialDiscretization ()-1;


		scalarVectorContainer_Type node_0( 2 );
		scalarVectorContainer_Type node_1( 2 );
		scalarVectorContainer_Type node_2( 2 );

		base_node t0(1);
		base_node t1(1);
		t0[ 0 ] = 0;
		t1[ 0 ] = 1;

		node_0[ 0 ].push_back( M_fractures [ 0 ]->getData ().getA ());
		node_0[ 1 ].push_back( M_fractures [ 0 ]->getData ().getB ());
		node_0[ 0 ].push_back( M_fractures [ 0 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_0[ 1 ].push_back( M_fractures [ 0 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		node_1[ 0 ].push_back( M_fractures [ 1 ]->getData ().getA ());
		node_1[ 1 ].push_back( M_fractures [ 1 ]->getData ().getB ());
		node_1[ 0 ].push_back( M_fractures [ 1 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_1[ 1 ].push_back( M_fractures [ 1 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		node_2[ 0 ].push_back( M_fractures [ 2 ]->getData ().getA ());
		node_2[ 1 ].push_back( M_fractures [ 2 ]->getData ().getB ());
		node_2[ 0 ].push_back( M_fractures [ 2 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_2[ 1 ].push_back( M_fractures [ 2 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		if( node_0 [ 0 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;

			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 0, 0 );
			M_intersection.setPoint ( 1, 0 );

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0. );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0. );
		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 0, nbDof0);
			M_intersection.setPoint ( 1, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
		}
		else if( node_0 [ 0 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 1 ]->getData ().imposeIntersection ( nbDof1 );
			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 1 ]->getData ().imposeIntersection ( nbDof1 );
			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}

		if( node_0 [ 0 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 0, 0);
			M_intersection.setPoint ( 2, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 0, nbDof0);
			M_intersection.setPoint ( 2, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 0 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 2 ]->getData ().imposeIntersection ( nbDof2 );
			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 2 ]->getData ().imposeIntersection ( nbDof2 );
			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
	}
	else if ( M_fractures.size() == 4 )
	{
		M_intersection.Resize(4);

		size_type nbDof0 = M_fractures [ 0 ]->getData().getSpatialDiscretization ()-1;
		size_type nbDof1 = M_fractures [ 1 ]->getData().getSpatialDiscretization ()-1;
		size_type nbDof2 = M_fractures [ 2 ]->getData().getSpatialDiscretization ()-1;
		size_type nbDof3 = M_fractures [ 3 ]->getData().getSpatialDiscretization ()-1;


		scalarVectorContainer_Type node_0( 2 );
		scalarVectorContainer_Type node_1( 2 );
		scalarVectorContainer_Type node_2( 2 );
		scalarVectorContainer_Type node_3( 2 );

		base_node t0(1);
		base_node t1(1);
		t0[ 0 ] = 0;
		t1[ 0 ] = 1;

		node_0[ 0 ].push_back( M_fractures [ 0 ]->getData ().getA ());
		node_0[ 1 ].push_back( M_fractures [ 0 ]->getData ().getB ());
		node_0[ 0 ].push_back( M_fractures [ 0 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_0[ 1 ].push_back( M_fractures [ 0 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		node_1[ 0 ].push_back( M_fractures [ 1 ]->getData ().getA ());
		node_1[ 1 ].push_back( M_fractures [ 1 ]->getData ().getB ());
		node_1[ 0 ].push_back( M_fractures [ 1 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_1[ 1 ].push_back( M_fractures [ 1 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		node_2[ 0 ].push_back( M_fractures [ 2 ]->getData ().getA ());
		node_2[ 1 ].push_back( M_fractures [ 2 ]->getData ().getB ());
		node_2[ 0 ].push_back( M_fractures [ 2 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_2[ 1 ].push_back( M_fractures [ 2 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		node_3[ 0 ].push_back( M_fractures [ 3 ]->getData ().getA ());
		node_3[ 1 ].push_back( M_fractures [ 3 ]->getData ().getB ());
		node_3[ 0 ].push_back( M_fractures [ 3 ]->getLevelSet ()->getData ( )->y_map( t0 ));
		node_3[ 1 ].push_back( M_fractures [ 3 ]->getLevelSet ()->getData ( )->y_map( t1 ));

		if( node_0 [ 0 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 0, 0);
			M_intersection.setPoint ( 1, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 1, 0);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 0 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( nbDof1 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 1 ]->getData ().imposeIntersection ( nbDof1 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 1 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}

		if( node_0 [ 0 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 2 ] = 0;

			M_intersection.setPoint ( 0, 0);
			M_intersection.setPoint ( 2, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 2, 0);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 0 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( nbDof2 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 2 ]->getData ().imposeIntersection ( nbDof2 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 2 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}

		if( node_0 [ 0 ][ 0 ] == node_3 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_3[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 3 ] = 0;
			M_intersection.setPoint ( 3, 0);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 3 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 3 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_3 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_3 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 3 ] = 0;
			M_intersection.setPoint ( 3, 0);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 3 ]->getData ().imposeIntersection ( 0 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 3 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );

		}
		else if( node_0 [ 0 ][ 0 ] == node_3 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_3 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 3 ] = nbDof3;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 3, nbDof3);
			M_intersection.setPoint ( 0, 0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( 0 );
			M_fractures [ 3 ]->getData ().imposeIntersection ( nbDof3 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 0, 0 );
			M_fractures [ 3 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}
		else if( node_0 [ 1 ][ 0 ] == node_3 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_3 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 3 ] = nbDof3;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 3, nbDof3);
			M_intersection.setPoint ( 0, nbDof0);

			M_fractures [ 0 ]->getData ().imposeIntersection ( nbDof0 );
			M_fractures [ 3 ]->getData ().imposeIntersection ( nbDof3 );

			M_fractures [ 0 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );
			M_fractures [ 3 ]->getData ().getFluxHandler( 0 )->update_UI ( 1, 0 );

		}

	}
	return;
}// updateNodes


void IntersectData::update_Ui ( const scalar_type& dt )
{
	/**
	 * u_I = u_I0 + 1/(|I|)*(F_0 + F_1 + F_2)
	 *
	 * dove u_I è il mio M_u0
	 */

	// devo calcolare tre flussi di Godunov, uno per ogni frattura
	scalarVector_Type Flux ( M_fractures.size() );

	scalar_type u0;

	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		if ( M_intersectionPoint [ i ] == 0 )
		{
			u0 = M_fractures [ i ]->getData().getFluxHandler ( 0 )->geta();
		}
		else
		{
			u0 = M_fractures [ i ]->getData().getFluxHandler ( 0 )->getb();
		}

		// calcolo il flusso di Godunov

		/***/

		size_type k = 0.;

		std::string monotone = M_fractures [ i ]->getData().getFluxHandler( k )->getMonotone();

		// calcolo il punto di massimo o di minimo, punto in cui si annulla la derivata
		scalar_type Us = 0.;

		size_type H = M_fractures [ i ]->getData().getFluxHandler( k )->getH();

		if( monotone == "false" )
		{
			Us = M_fractures [ i ]->getData().getFluxHandler( k )->getUs();
		}

		if ( monotone == "false" )
		{
			scalar_type Sl;
			scalar_type Sr;


			if( H == 2.0 )		// caso A
			{
				if ( M_intersectionPoint [ i ] == 0 )
				{
					Sl = std::max( M_u0, Us );
					Sr = std::min( u0, Us );
				}
				else
				{
					Sl = std::max( u0, Us );
					Sr = std::min( M_u0, Us );
				}


				scalar_type F_ul = M_fractures [ i ]->getData().feval_scal( Sl, k, M_intersectionPoint [ i ]-1 );
				scalar_type F_ur = M_fractures [ i ]->getData().feval_scal( Sr, k, M_intersectionPoint [ i ]-1 );

				Flux[ i ] = std::max( F_ul, F_ur );
			}

			else		// caso B
			{
				if ( M_intersectionPoint [ i ] == 0 )
				{
					Sl = std::min( M_u0, Us );
					Sr = std::max( u0, Us );
				}
				else
				{
					Sl = std::min( u0, Us );
					Sr = std::max( M_u0, Us );
				}


				scalar_type F_ul = M_fractures [ i ]->getData().feval_scal( Sl, k, M_intersectionPoint [ i ] );
				scalar_type F_ur = M_fractures [ i ]->getData().feval_scal( Sr, k, M_intersectionPoint [ i ] );

				Flux[ i ] = std::min( F_ul, F_ur );
			}
		}
		else
		{
			scalar_type F_ul = M_fractures [ i ] ->getData().feval_scal( u0, k, M_intersectionPoint [ i ] );
			scalar_type F_ur = M_fractures [ i ] ->getData().feval_scal( M_u0, k, M_intersectionPoint [ i ] );

			scalar_type Sl = u0;
			scalar_type Sr = M_u0;

			if ( M_intersectionPoint [ i ] == 0 )
			{
				F_ur = M_fractures [ i ] ->getData().feval_scal( u0, k, M_intersectionPoint [ i ] );
				F_ul = M_fractures [ i ] ->getData().feval_scal( M_u0, k , M_intersectionPoint [ i ]);

				Sl = M_u0;
				Sr = u0;
			}

			scalar_type s = ( F_ur - F_ul )*( Sr - Sl );

			if ( s > 0 )
			{
				Flux[ i ] = F_ul;
			}
			else
			{
				Flux[ i ] = F_ur;
			}

		}
	}

					// ***********  //
					// modifico qui //
					// ***********  //

	// a questo punto calcolo la nuova U_I
	scalar_type U0_old = M_u0;

	if ( M_intersectionPoint [0] == 0 )
	{
		//if ( ( M_fractures [ 0 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 0 ] =  - M_fractures[ 0 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 0 ] =  - M_fractures[ 0 ]->getData().getSI()*M_fractures[0]->getData().getThickness ();
		}*/
	}
	else
	{
		//if (( M_fractures [ 0 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 0 ] = M_fractures[ 0 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 0 ] =  M_fractures[ 0 ]->getData().getSI()*M_fractures[0]->getData().getThickness ();
		}*/
	}

	if ( M_intersectionPoint [1] == 0 )
	{
		//if ( (M_fractures [ 1 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 1 ] =  - M_fractures[ 1 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 1 ] =  - M_fractures[ 1 ]->getData().getSI()*M_fractures[1]->getData().getThickness ();
		}*/
	}
	else
	{
		//if ( (M_fractures [ 1 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 1 ] =  M_fractures[ 1 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 1 ] = M_fractures[ 1 ]->getData().getSI()*M_fractures[1]->getData().getThickness ();
		}*/
	}

	if ( M_intersectionPoint [2] == 0 )
	{
		//if ( (M_fractures [ 2 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 2 ] =  - M_fractures[ 2 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 2 ] =  - M_fractures[ 2 ]->getData().getSI()*M_fractures[2]->getData().getThickness ();
		}*/
	}
	else
	{

		//if ( (M_fractures [ 2 ] -> getData().getFluxHandler( 0 ))->getU() != "v")
		{
			Flux [ 2 ] =  M_fractures[ 2 ]->getData().getSI();
		}/*
		else
		{
			Flux [ 2 ] =  M_fractures[ 2 ]->getData().getSI()*M_fractures[2]->getData().getThickness ();
		}*/

	}

	M_u0 = U0_old - ( dt )/( M_measure ) *( - Flux [ 0 ] - Flux [ 1 ] - Flux [ 2 ] );

	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		M_fractures [ i ] -> getData().getFluxHandler( 0 )->update_UI ( M_intersectionPoint [ i ], M_u0 );
	}

	return;

}// update_Ui



IntersectData & IntersectData::operator =(const IntersectData& t)
{
	M_intersectionPoint = t.M_intersectionPoint;
	M_fractures = t.M_fractures;
	M_u0 = t.M_u0;
	M_measure = t.M_measure;

	M_intersection = t.M_intersection;

	return *this;

}// operator =
