
#include "../include/Geometry.h"

/**************************************************************************/
/*  Geometry.h															  */
/*  Libreria che costruisce il triangolo 2d dell'intersezione             */
/**************************************************************************/


FractureEnd::FractureEnd(PointData& a , scalar_type t)
{
	M_endPoint = a;

	M_thickness = t;
}//costruttore FractureEnd
 
Intersection::Intersection( FractureEndContainer_Type const & fractureEnd, PointData const & intersectionPoint) : 
								M_intersection (intersectionPoint)
{
	M_fractures.clear();

	for (size_type i = 0; i < fractureEnd.size(); i++) 
	{
		M_fractures.push_back( fractureEnd [ i ]);
	}
	
	Vector2d ttmp;
	Vector2d ntmp;
		
	if (fractureEnd.size() == 3)
	{
	
		for ( size_type i = 0; i < 3; ++i )
		{
			// get tangent at the end
			ttmp(0) = intersectionPoint.x()-M_fractures[i].getPoint ().x();

			ttmp(1) = intersectionPoint.y()-M_fractures[i].getPoint ().y();

			ttmp.normalize();
	
			M_tangents [ i ] = ttmp;
			ntmp(0) = -ttmp(1); 
			ntmp(1) = ttmp(0); 
	
			M_normals [ i ] = ntmp;
		}
		
		M_intersectionTriangle = this->computeIntersectionTriangle();
	}
	else if ( fractureEnd.size() == 4 )
	{
		for ( size_type i = 0; i < 4; ++i )
		{
			// get tangent at the end
			ttmp(0) = intersectionPoint.x()-M_fractures[i].getPoint ().x();

			ttmp(1) = intersectionPoint.y()-M_fractures[i].getPoint ().y();

			ttmp.normalize();

			M_tangents [ i ] = ttmp;
			ntmp(0) = -ttmp(1);
			ntmp(1) = ttmp(0);

			M_normals [ i ] = ntmp;
		}

		M_intersectionQuadrilater = this->computeIntersectionQuadrilater();
	}
	
}// costruttore intersezione


void Intersection::setIntersection( FracturePtrContainer_Type& M_FracturesSet )
{
	if ( M_FracturesSet.size() == 3 )
	{
		setTriangleIntersection( M_FracturesSet );

	} 
	else if ( M_FracturesSet.size() == 4 )
	{
		setQuadrilaterIntersection( M_FracturesSet );
	}

	return; 
	
}// setIntersection

void Intersection::setQuadrilaterIntersection( FracturePtrContainer_Type& M_FracturesSet ) 
{

	base_node node0(2);
	base_node node1(2);
	base_node node2(2);
	base_node node3(2);
	
	base_node nodeI(2);
	
	//trucchetto xk y_map vuole i base_node
	base_node tmp0(1);
	base_node tmp1(1);
	tmp0[0] = 0.;
	tmp1[0] = 1.;
	
	FractureHandlerPtr_Type f0, f1, f2, f3;

	f0 = M_FracturesSet [ 0 ];
	f1 = M_FracturesSet [ 1 ];
	f2 = M_FracturesSet [ 2 ];
	f3 = M_FracturesSet [ 3 ];

	node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp0);
	node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp0);
	node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp0);
	node3[0] = M_FracturesSet[3]->getLevelSet()->getData()->x_map(tmp0);
	nodeI[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);


	node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp0);
	node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp0);
	node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp0);
	node3[1] = M_FracturesSet[3]->getLevelSet()->getData()->y_map(tmp0);
	nodeI[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);


	if ( node0[0] == node1[0] && node0[1] == node1[1] )
	{
		nodeI[0] = node0[0];
		nodeI[1] = node0[1];

		node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
		node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);

		node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
		node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);

		if ( node2[0] == nodeI[0] && node2[1] == nodeI[1] )
		{
			node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

			node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
		}

		if ( node3[0] == nodeI[0] && node3[1] == nodeI[1] )
		{
			node3[0] = M_FracturesSet[3]->getLevelSet()->getData()->x_map(tmp1);

			node3[1] = M_FracturesSet[3]->getLevelSet()->getData()->y_map(tmp1);
		}
	}
	else if ( node0[0] == node2[0] && node0[1] == node2[1] )
	{
		nodeI[0] = node0[0];
		nodeI[1] = node0[1];

		node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
		node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

		node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
		node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);

		if ( node1[0] == nodeI[0] && node1[1] == nodeI[1] )
		{
			node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);

			node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
		}

		if ( node3[0] == nodeI[0] && node3[1] == nodeI[1] )
		{
			node3[0] = M_FracturesSet[3]->getLevelSet()->getData()->x_map(tmp1);

			node3[1] = M_FracturesSet[3]->getLevelSet()->getData()->y_map(tmp1);
		}
	}
	else if ( node0[0] == node3[0] && node0[1] == node3[1] )
	{
		nodeI[0] = node0[0];
		nodeI[1] = node0[1];

		node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
		node3[0] = M_FracturesSet[3]->getLevelSet()->getData()->x_map(tmp1);

		node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
		node3[1] = M_FracturesSet[3]->getLevelSet()->getData()->y_map(tmp1);

		if ( node1[0] == nodeI[0] && node1[1] == nodeI[1] )
		{
			node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);

			node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
		}

		if ( node2[0] == nodeI[0] && node2[1] == nodeI[1] )
		{
			node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

			node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
		}
	}
	else if ( node1[0] == nodeI[0] && node1[1] == nodeI[1] )
	{
		node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);

		node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
	}
	//else
	if ( node2[0] == nodeI[0] && node2[1] == nodeI[1] )
	{
		node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

		node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
	}
	//else
	if ( node3[0] == nodeI[0] && node3[1] == nodeI[1] )
	{
		node3[0] = M_FracturesSet[3]->getLevelSet()->getData()->x_map(tmp1);

		node3[1] = M_FracturesSet[3]->getLevelSet()->getData()->y_map(tmp1);
	}


	//Estremi liberi delle fratture +  punto di intersezione
	PointData p0(node0[0], node0[1]);
	PointData p1(node1[0], node1[1]);
	PointData p2(node2[0], node2[1]);
	PointData p3(node3[0], node3[1]);
	
	PointData pi(nodeI[0], nodeI[1]);
	
	//Spessori delle mie fratture
	scalar_type t0, t1, t2, t3;
	t0 = M_FracturesSet[0]->getData().getThickness();
	t1 = M_FracturesSet[1]->getData().getThickness();
	t2 = M_FracturesSet[2]->getData().getThickness();
	t3 = M_FracturesSet[3]->getData().getThickness();
	
	FractureEnd fracture0(p0, t0);
	FractureEnd fracture1(p1, t1);
	FractureEnd fracture2(p2, t2);
	FractureEnd fracture3(p3, t3);

	
	FractureEndContainer_Type fractureEnd;
	fractureEnd.push_back(fracture0);
	fractureEnd.push_back(fracture1);
	fractureEnd.push_back(fracture2);
	fractureEnd.push_back(fracture3);

	
	Intersection tmp(fractureEnd, pi);


	this->M_fractures = tmp.M_fractures;
	this->M_intersection = tmp.M_intersection;
	
	this->M_tangents[0] = tmp.M_tangents[0];
	this->M_tangents[1] = tmp.M_tangents[1];
	this->M_tangents[2] = tmp.M_tangents[2];
	this->M_tangents[3] = tmp.M_tangents[3];
	
	this->M_normals[0] = tmp.M_normals[0];
	this->M_normals[1] = tmp.M_normals[1];
	this->M_normals[2] = tmp.M_normals[2];
	this->M_normals[3] = tmp.M_normals[3];
	
	this->M_intersectionQuadrilater = tmp.M_intersectionQuadrilater;
	
	std::cout << tmp.M_intersectionQuadrilater << std::endl;


	return;

}// setQuadrilaterIntersection


// posso modificarla
void Intersection::setTriangleIntersection(	FracturePtrContainer_Type& M_FracturesSet ) 
{
	base_node node0(2);
	base_node node1(2);
	base_node node2(2);
	base_node nodeI(2);
	
	base_node tmp0(1);
	base_node tmp1(1);
	tmp0[0] = 0.;
	tmp1[0] = 1.;


	node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp0);
	node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp0);
	node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp0);
	nodeI[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);


	node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp0);
	node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp0);
	node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp0);
	nodeI[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);

	if ( node0[0] == node1[0] && node0[1] == node1[1] )
	{
		if ( node2[0] == node0[0] && node2[1] == node0[1] )
		{
			nodeI[0] = node0[0];
			nodeI[1] = node0[1];

			node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
			node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);
			node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

			node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
			node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
			node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
		}
		else
		{
			nodeI[0] = node0[0];
			nodeI[1] = node0[1];

			node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
			node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);

			node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
			node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
		}

	}
	else if ( node2[0] == node0[0] && node2[1] == node0[1] )
	{
		nodeI[0] = node0[0];
		nodeI[1] = node0[1];

		node0[0] = M_FracturesSet[0]->getLevelSet()->getData()->x_map(tmp1);
		node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

		node0[1] = M_FracturesSet[0]->getLevelSet()->getData()->y_map(tmp1);
		node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
	}
	else if ( node2[0] == node1[0] && node2[1] == node1[1] )
	{
		node1[0] = M_FracturesSet[1]->getLevelSet()->getData()->x_map(tmp1);
		node2[0] = M_FracturesSet[2]->getLevelSet()->getData()->x_map(tmp1);

		node1[1] = M_FracturesSet[1]->getLevelSet()->getData()->y_map(tmp1);
		node2[1] = M_FracturesSet[2]->getLevelSet()->getData()->y_map(tmp1);
	}


	//Estremi liberi delle fratture +  punto di intersezione
	PointData p0(node0[0], node0[1]);
	PointData p1(node1[0], node1[1]);
	PointData p2(node2[0], node2[1]);
	PointData pi(nodeI[0], nodeI[1]);
	

	//Spessori delle mie fratture
	scalar_type t0, t1, t2;
	t0 = M_FracturesSet[0]->getData().getThickness();
	t1 = M_FracturesSet[1]->getData().getThickness();
	t2 = M_FracturesSet[2]->getData().getThickness();
	
	FractureEnd f0(p0, t0);
	FractureEnd f1(p1, t1);
	FractureEnd f2(p2, t2);
	
	FractureEndContainer_Type fractureEnd;
	fractureEnd.push_back(f0);
	fractureEnd.push_back(f1);
	fractureEnd.push_back(f2);
	
	Intersection tmp(fractureEnd, pi);
	

	this->M_fractures = tmp.M_fractures;
	this->M_intersection = tmp.M_intersection;
	
	this->M_tangents[0] = tmp.M_tangents[0];
	this->M_tangents[1] = tmp.M_tangents[1];
	this->M_tangents[2] = tmp.M_tangents[2];
	
	this->M_normals[0] = tmp.M_normals[0];
	this->M_normals[1] = tmp.M_normals[1];
	this->M_normals[2] = tmp.M_normals[2];
	
	this->M_intersectionTriangle = tmp.M_intersectionTriangle;
	
	std::cout << std::endl;
	std::cout << tmp.M_intersectionTriangle << std::endl;

	return;

}// setTriangleIntersection


TriangleData const & Intersection::computeIntersectionTriangle()
{
	PointData point;
	const scalar_type tol=1.e-5;
	scalar_type s(0.);


	for ( size_type i = 0; i < 3; ++i )
	{
		size_type j = (i +1 ) % 3;
		scalar_type ninj = M_normals[i].dot(M_normals[j] );
		scalar_type nitj = M_normals[i].dot(M_tangents[j]);
		scalar_type titj = M_tangents[i].dot(M_tangents[j]);

		if(std::fabs(nitj)<tol)
		{
			point = M_intersection +
					PointData ( M_normals [ i ] ( 0 ), M_normals [ i ] ( 1 ) ) * ( 0.5 * ( M_fractures [ j ].getThickness() + M_fractures [ i ].getThickness() ));
		}
		else
		{
			// parametric coordinate
			s = 0.5 * ( M_fractures [ i ].getThickness() + M_fractures [ j ].getThickness() * ninj ) / nitj;

			// the ith point is Pji in the note
			point = M_intersection + ( PointData ( M_tangents [ j ] ( 0 ), M_tangents [ j ] ( 1 ) ) * s )-
									( PointData ( M_normals [ j ] ( 0 ), M_normals [ j ] ( 1 ) ) * ( 0.5 * M_fractures [ j ].getThickness() ) );
		}

		M_intersectionTriangle.setPoint(i,point);

	}

	return M_intersectionTriangle;

}// computeIntersectionTriangle

QuadrilaterData const & Intersection::computeIntersectionQuadrilater()
{
	PointData point;
	const scalar_type tol=1.e-5;
	scalar_type s(0.);
	
	for (size_type i=0; i<4;++i)
	{
		size_type j = (i +1 )%4;
		size_type k, p;
		
		if( i == 0 || i ==2 )
		{
			k = 0;
			p = 1;
		}
		else
		{
			k = 1;
			p = 0;
		}
		scalar_type ninj = M_normals[i].dot(M_normals[j] );
		scalar_type nitj = M_normals[i].dot(M_tangents[j]);
		 
		if(std::fabs(nitj)<tol)
		{
			point = M_intersection + 
					PointData ( M_normals [ i ] ( 0 ), M_normals [ i ] ( 1 ) ) * ( 0.5 * ( M_fractures [ p ].getThickness() + M_fractures [ k ].getThickness() ));
		}
		else
		{
			// parametric coordinate
			s    = 0.5 * ( M_fractures [ p ].getThickness() * ninj + M_fractures [ k ].getThickness() )/nitj;
			
			// the ith point is Pji in the note
			point = M_intersection + ( PointData ( M_tangents [ j ] ( 0 ), M_tangents [ j ] ( 1 ) ) * s )- 
									( PointData ( M_normals [ j ] ( 0 ), M_normals [ j ] ( 1 ) ) * ( 0.5 * M_fractures [ p ].getThickness() ) );
		}
	
		M_intersectionQuadrilater.setPoint(i,point);
		
	}
	
	return M_intersectionQuadrilater;

}// computeIntersectionQuadrilater




void Intersection::Resize ( const size_type& i )
{
	M_intersectionPoint.clear();
	M_intersectionPoint.resize(i);

	return;
}// Resize


void Intersection::setPoint ( const size_type& i, const size_type& j)
{
	M_intersectionPoint [ i ] = j;

	return;
}// setPoint


void Intersection::sortFractures ( FracturePtrContainer_Type& M_FracturesSet )
{
	FractureHandlerPtr_Type f0 = M_FracturesSet [ 0 ];
	FractureHandlerPtr_Type f1 = M_FracturesSet [ 1 ];
	FractureHandlerPtr_Type f2 = M_FracturesSet [ 2 ];

	if ( isPos ( f0, f1) )
	{
		if ( isPos ( f0, f2 ) )
		{
			if ( isPos ( f1, f2 ) )
			{
				return;
			}
			else
			{
				// scambio 1 e 2
				FractureHandlerPtr_Type tmp = f1;
				M_FracturesSet [ 1 ] = M_FracturesSet [ 2 ];
				M_FracturesSet [ 2 ] = tmp;
				return;
			}
		}
		else
		{
			return;
		}
	}
	else	// sistemo da qui
	{
		if ( isPos ( f0, f2 ) )
		{
			// scambio 1 e 2
			FractureHandlerPtr_Type tmp = f1;

			M_FracturesSet [ 1 ] = M_FracturesSet [ 2 ];
			M_FracturesSet [ 2 ] = tmp;
			return;
		}
		else if ( !isPos( f1, f2 ) )
		{
			// scambio 1 e 2
			FractureHandlerPtr_Type tmp = f1;

			M_FracturesSet [ 1 ] = M_FracturesSet [ 2 ];
			M_FracturesSet [ 2 ] = tmp;
			return;
		}
		else
		{
			return;
		}
	}

	return;
}// sortFractures



bool Intersection::isPos ( FractureHandlerPtr_Type& fracture, FractureHandlerPtr_Type& otherfracture)
{
	base_node n0(2);

	n0 [ 0 ] = 0;
	n0 [ 1 ] = otherfracture->getLevelSet()->getData()->y_map ( n0 );
	n0 [ 0 ] = otherfracture->getLevelSet()->getData()->x_map ( n0 );

	base_node n1(2);

	n1 [ 0 ] = 1;
	n1 [ 1 ] = otherfracture->getLevelSet()->getData()->y_map ( n1 );
	n1 [ 0 ] = otherfracture->getLevelSet()->getData()->x_map ( n1 );

	if ( fracture->getLevelSet()->getData()->ylevelSetFunction ( n0 ) >= 0 && fracture->getLevelSet()->getData()->levelSetFunction ( n1 ) >= 0 )
	{
		return true;
	}
	else 
	{
		return false;
	}
	
}// isPos
