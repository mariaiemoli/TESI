
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

	scalar_type s = 0.;
	scalar_type t = 0.;

	FracturePtrContainer_Type fractures_copy = fractures;

	M_intersection.setIntersection ( fractures_copy );

	for ( size_type i = 0; i < fractures.size(); i++ )
	{
		t = fractures [ i ]->getMeshFEM().point_of_basic_dof( M_intersectionPoint [ i ] )[ 0 ];
		//s = s + 1./( fractures [ i ]->getData ().getCi ( t )  );
		s = s + ( fractures [ i ]->getData ().getCi ( t )  );
	}

	// per la condizione iniziale nel punto di intersezione prendo la media armonica del valore della ci
	// nei punti più prossimi all'intersezione delle fratture

	//M_u0 = fractures.size()/s;
	M_u0 = s/fractures.size();

	M_u0 = 1.;

	// aggiorno il dato al bordo per le fratture che si intersecano
	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		update_Bc ( i );
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
			M_intersection.setPoint ( 0, 0);
			M_intersection.setPoint ( 1, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 0, nbDof0);
			M_intersection.setPoint ( 1, 0);
		}
		else if( node_0 [ 0 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, 0);

		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, nbDof0);

		}

		if( node_0 [ 0 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 0, 0);
			M_intersection.setPoint ( 2, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 0, nbDof0);
			M_intersection.setPoint ( 2, 0);
		}
		else if( node_0 [ 0 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, nbDof0);
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
		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 1 ] = 0;
			M_intersection.setPoint ( 1, 0);
			M_intersection.setPoint ( 0, nbDof0);
		}
		else if( node_0 [ 0 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_1 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_1 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 1 ] = nbDof1;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 1, nbDof1);
			M_intersection.setPoint ( 0, nbDof0);
		}

		if( node_0 [ 0 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 2 ] = 0;
		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 2 ] = 0;
			M_intersection.setPoint ( 2, 0);
			M_intersection.setPoint ( 0, nbDof0);
		}
		else if( node_0 [ 0 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_2 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_2 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 2 ] = nbDof2;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 2, nbDof2);
			M_intersection.setPoint ( 0, nbDof0);
		}

		if( node_0 [ 0 ][ 0 ] == node_3 [ 0 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_3[ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = 0;
			M_intersectionPoint [ 3 ] = 0;
			M_intersection.setPoint ( 3, 0);
			M_intersection.setPoint ( 0, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_3 [ 0 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_3 [ 0 ][ 1 ] )
		{
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersectionPoint [ 3 ] = 0;
			M_intersection.setPoint ( 3, 0);
			M_intersection.setPoint ( 0, nbDof0);

		}
		else if( node_0 [ 0 ][ 0 ] == node_3 [ 1 ][ 0 ] && node_0 [ 0 ][ 1 ] == node_3 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 3 ] = nbDof3;
			M_intersectionPoint [ 0 ] = 0;
			M_intersection.setPoint ( 3, nbDof3);
			M_intersection.setPoint ( 0, 0);
		}
		else if( node_0 [ 1 ][ 0 ] == node_3 [ 1 ][ 0 ] && node_0 [ 1 ][ 1 ] == node_3 [ 1 ][ 1 ] )
		{
			M_intersectionPoint [ 3 ] = nbDof3;
			M_intersectionPoint [ 0 ] = nbDof0;
			M_intersection.setPoint ( 3, nbDof3);
			M_intersection.setPoint ( 0, nbDof0);
		}

	}
	return;
}// updateNodes


void IntersectData::update_Ui ( const scalar_type& landa )
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
			u0 = M_fractures [ i ]->getData().getFluxHandler ( 0 )->getUin();
		}
		else
		{
			u0 = M_fractures [ i ]->getData().getFluxHandler ( 0 )->getUout();
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
			if( H == 2.0 )		// caso A
			{
				scalar_type m1 = std::max( u0, Us );
				scalar_type m2 = std::min( M_u0, Us );

				scalar_type F_ul = M_fractures [ i ]->getData().feval_scal( m1, k );
				scalar_type F_ur = M_fractures [ i ]->getData().feval_scal( m2, k );

				std::cout << " frattura  " << i << " F_ul  " << F_ul << "   F_ur  " << F_ur << std::endl;
				Flux[ i ] = std::max( F_ul, F_ur );
			}

			else		// caso B
			{
				scalar_type m1 = std::min( u0, Us );
				scalar_type m2 = std::max( M_u0, Us );

				scalar_type F_ul = M_fractures [ i ] ->getData().feval_scal( m1, k );
				scalar_type F_ur = M_fractures [ i ] ->getData().feval_scal( m2, k );

				Flux[ i ] = std::min( F_ul, F_ur );
			}
		}
		else
		{
			scalar_type F_ul = M_fractures [ i ] ->getData().feval_scal( u0, k );
			scalar_type F_ur = M_fractures [ i ] ->getData().feval_scal( M_u0, k );

			scalar_type s = (F_ur - F_ul )*(M_u0 - u0 );


			if ( i == 1 ) //( s >= 0)
			{
				Flux[ i ] = F_ul;
			}
			else
			{
				Flux[ i ] = F_ur;
			}

		}
		/***/
	}

	// a questo punto calcolo la nuova U_I
	scalar_type U0_old = M_u0;
	/*
	size_type H0 = M_fractures[ 0 ]->getData().getFluxHandler( 0 )->getH();
	size_type H1 = M_fractures[ 1 ]->getData().getFluxHandler( 0 )->getH();
	size_type H2 = M_fractures[ 2 ]->getData().getFluxHandler( 0 )->getH();
	*/


	scalar_type H0 = M_fractures[ 0 ]->getData().getThickness ();
	scalar_type H1 = M_fractures[ 1 ]->getData().getThickness ();
	scalar_type H2 = M_fractures[ 2 ]->getData().getThickness ();

	//M_measure = - M_measure;

	/*
	std::cout << " Flux [ 0 ] " << Flux [ 0 ] << std::endl;
	std::cout << " Flux [ 1 ] " << Flux [ 1 ] << std::endl;
	std::cout << " Flux [ 2 ] " << Flux [ 0 ] << std::endl;
	*/

	//std::cout << "  ( landa )/( M_measure ) *(  - Flux [ 0 ] - Flux [ 1 ] - Flux [ 2 ] ) " << ( landa )/( M_measure ) *(  - Flux [ 0 ] - Flux [ 1 ] - Flux [ 2 ] ) << std::endl;

	M_u0 = U0_old - ( landa )/( M_measure ) *(  - Flux [ 0 ] - Flux [ 1 ] - Flux [ 2 ] );

	//M_u0 = U0_old - ( landa )/( M_measure ) *( - Flux [ 0 ]*H0 + Flux [ 1 ]*H1 - Flux [ 2 ]*H2 );

	//std::cout << " M_u0: " << M_u0 << std::endl;

	//M_u0 = U0_old - ( landa )/( M_measure ) *( Flux [ 0 ]*H0 + H1*Flux [ 1 ] + H2*Flux [ 2 ] );
	//M_u0 = U0_old - ( landa )*( Flux [ 0 ] + Flux [ 1 ] + Flux [ 2 ] );

	for ( size_type i = 0; i < M_fractures.size(); i++ )
	{
		update_Bc ( i );
	}

	return;

}// update_Ui



void IntersectData::update_Bc ( const size_type& i )
{
	size_type nbDof = M_fractures [ i ]->getData ().getSpatialDiscretization ()-1;

	if ( M_intersectionPoint [ i ] == nbDof  )
	{
		M_fractures [ i ] -> getData().update_Bc ( nbDof -1, M_u0);
	}
	else
	{
		M_fractures [ i ] -> getData().update_Bc ( 0, M_u0);
	}

	return;
}// update_Bc




IntersectData & IntersectData::operator =(const IntersectData& t)
{
	M_intersectionPoint = t.M_intersectionPoint;
	M_fractures = t.M_fractures;
	M_u0 = t.M_u0;
	M_measure = t.M_measure;

	M_intersection = t.M_intersection;

	return *this;

}// operator =
