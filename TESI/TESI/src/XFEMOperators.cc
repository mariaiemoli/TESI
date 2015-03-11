
#include "../include/XFEMOperators.h"

/**************************************************************************/
/*  XFEMOperators.cc												  	  */
/*  Libreria che definisce le forme lineari e bilineari per il problema	  */
/*  di Darcy             												  */
/**************************************************************************/

namespace getfem
{

// Defining unit normal on a level set ------------------------------------

level_set_unit_normal::level_set_unit_normal ( const getfem::mesh_fem& mf_, const scalarVector_Type& U_ ) :
    mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()), gradU(1, N)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = short_type(N);
    mf.extend_vector(U_, U);
}

void level_set_unit_normal::compute ( getfem::fem_interpolation_context& ctx, bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    
    gmm::copy( gmm::sub_vector( U, gmm::sub_index( mf.ind_basic_dof_of_element( cv ) ) ), coeff );

    ctx.pf()->interpolation_grad( ctx, coeff, gradU, 1);
    
    scalar_type norm = gmm::vect_norm2( gmm::mat_row( gradU, 0 ));
    
    for ( size_type i = 0; i < N; ++i )
    {
        t [ i ] = gradU(0, i) / norm;
    }
    
    return;
}



//matrice A11 per la frattura non intersecata
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalar_type& gammaU,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const sizeVector_Type& ExtBoundary,
                  const size_type& uncutRegionFlag )
{
    const size_type shiftVelocity = fracture->getMeshFEM2().nb_dof();
    const size_type shiftData = fracture->getMeshFEM().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize ( M_, shiftVelocity, shiftVelocity );
    gmm::clear ( M_ );

    scalarVector_Type etaGammaUinvh ( shiftData );

    for ( size_type i = 0; i < shiftData; ++i )
    {
        etaGammaUinvh [ i ] = invKTangentialInterpolated [ i ] * gammaU * fracture->getInverseMeshSize(i);
    }

    // Volume integration
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF ( shiftMapFactor, 0. );

    generic_assembly assem;

    /**
     * getfem::generic_assembly assem;
     *
     * assem.push_im(mim);
     * assem.push_mf(mf);
     * assem.push_mf(mfdata);
     * assem.push_data(F);
     * assem.push_vec(B);
     *
     * assem.set("Z=data(#2);"
     * 			 "V(#1)+=comp(Base(#1).Base(#2))(:,j).Z(j);");
     *
     * assem.assembly();
     *
     * The first instructions declare the object, and set the data that it will use: a mesh_im object which holds the integration
     * methods, two mesh_fem objects, the input data F, and the destination vector B.
     *
     * The input data is the vector F , defined on mfdata.
     * One wants to evaluate sum(j){ f_j* int_Ω (φ_i * ψ_j). The instruction must be seen as something that will be executed for each convex cv
     * of the mesh.
     * The terms #1 and #2 refer to the first mesh_fem and the second one (i.e. mf and mfdata).
     * The instruction Z=data(#2); means that for each convex, the “tensor” Z will receive the values of the first data argument provided
     * with push_data, at indexes corresponding to the degrees of freedom attached to the convex of the second (#2) mesh_fem
     * (here, Z = F[mfdata.ind_dof_of_element(cv)].
     * The part V(#1)+=... means that the result of the next expression will be accumulated into the output vector (provided with push_vec).
     * Here again, #1 means that we will write the result at indexes corresponding to the degrees of freedom of the current convex with
     * respect to the first (#1) mesh_fem.
     *
     * The right hand side comp(Base(#1).Base(#2))(:,j).Z(j) contains two operations.
     * The first one is a computation of a tensor on the convex: comp(Base(#1).Base(#2)) is evaluated as a 2-dimensions tensor, int(φ_i*ψ_j) ,
     * for all degrees of freedom i of mf and j of mfdata attached to the current convex.
     * The next part is a reduction operation, C(:,j).Z(j): each named index (here j) is summed, i.e. the result is sum(j){ c_(i,j)*z_j }.
     *
     * The integration method used inside comp(Base(#1).Base(#2)) is taken from mim.
     *
     */

    if ( fracture->getMeshFEM2().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        /*
         *  definisce la forma bilineare:
         *
         *  		a_i(u,w) = (eta_i * u, w)_L2
         *
         *  u velocità
         *  w funzione test per la velocità
         *
         *  #1 velocità
         *  #2 pressione
         *
         */
        assem.set("w=data$1(#2);" "q=data$2(#2);"
        		  "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
        		  "M(#1,#1)+=a(:, :, i,k).w(i).q(k);");
    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / (fracture->getMagnificationMapFactor1(i)
                    * fracture->getMagnificationMapFactor2(i));
        }
        assem.set("w=data(#2);" "q=data$2(#2);"
            "a=comp(vBase(#1).vBase(#1).Base(#2).Base(#2));"
            "M(#1,#1)+=a(:,i,:,i,j,k).w(j).q(k);");
    }

    // Assegno il metodo di integrazione su M_mediumMesh
    assem.push_mi(fracture->getIntegrationMethod2());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh
    assem.push_mf(fracture->getMeshFEM2());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh per i coefficienti
    assem.push_mf(fracture->getMeshFEM());
    assem.push_mf(fracture->getMeshFEM());

    // Assegno i coefficienti
    assem.push_data(invKTangentialInterpolated);
    assem.push_data(invF);

    // Definisco la matrice dove salare i risultati
    assem.push_mat(M_);

    // Calcolo la matrice
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator a(volume)      [OK]" << endl;

    // Boundary integration for the fracture
    gmm::clear(M_);

    getfem::generic_assembly assem_surf;

    /*
     * tratta il termine di bordo non soggetto a condizione al contorno di dirichlet per la pressione
     *
     */
    assem_surf.set("gamma=data$1(#2);"
    			   "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
    			   "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assegno il metodo di integrazione su M_mediumMesh
    assem_surf.push_mi(fracture->getIntegrationMethod2());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh
    assem_surf.push_mf(fracture->getMeshFEM2());

    // Assegno lo spazio degli elementi finiti su M_mediumMesh per i coefficienti
    assem_surf.push_mf(fracture->getMeshFEM());

    // Assegno i coefficienti
    assem_surf.push_data(etaGammaUinvh);

    // Definisco la matrice dove salare i risultati
    assem_surf.push_mat(M_);

    // Assemblo la matrice su ogni sottoregione
    for ( size_type bndID = 0; bndID < ExtBoundary.size(); bndID++ )
    {
        assem_surf.assembly(
                fracture->getMeshFEM2().linked_mesh().get_mpi_sub_region(
                        ExtBoundary [ bndID ]));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) += M_(i, j);
        }
    }
    cout << "DARCY :: operator a(surface)     [OK]" << endl;

    return;

} // darcy_A11F


//matrice A12 per la frattura non intersecata
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& uncutRegionFlag )
{
    const size_type velocityShift = fracture->getMeshFEM2().nb_dof();
    const size_type pressureShift = fracture->getMeshFEM().nb_dof();

    sparseMatrix_Type M_;

    gmm::resize(M_, velocityShift, pressureShift);
    gmm::clear(M_);

    // Volume integration
    getfem::generic_assembly assem;

    if ( fracture->getMeshFEM2().get_qdim() == 1 )
    {
    	/*
    	 * definisce la forma bilineare
    	 *
    	 * 		b_i(q,w) = - (q, div(w))_L2
    	 *
    	 * 	w funzione test per la velocità
    	 * 	q funzione test per la pressione
    	 *
    	 */
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:, i,i,:);");
    }
    else
    {
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:,i,i, :);");
    }

    // Assign the M_mediumMesh integration method
    assem.push_mi(fracture->getIntegrationMethod2());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(fracture->getMeshFEM2());
    assem.push_mf(fracture->getMeshFEM());

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < velocityShift; ++i )
    {
        for ( size_type j = 0; j < pressureShift; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator b(volumic)     [OK]" << endl;

    return;

} // darcy_A12F



void setAup_i ( sparseMatrixPtr_Type& Aup_i,
				size_type id, size_type id_i, size_type id_j, size_type id_k,
				scalarVector_Type& DOF, scalarVector_Type& DOF_v,
				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,
				Matrix4d T, const size_type Index,
				scalar_type s )
{
	size_type id0, id1, id2;

	id0 = id_i;
	id1 = id_j;
	id2 = id_k;

	orderId( id0, id1, id2 );


	if( id != 3 )
	{
		(*Aup_i) ( 0 , shiftIntersect[ id_i ] + DOF_v[ id ] )  = 1.;
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( id , 0 );
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( id , 1 );
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = 1.*T( id , 2 );
		(*Aup_i) ( 0 , Index ) = -1.*( T( id , 0 ) + T( id , 1 ) + T( id , 2 ) );
	}
	else
	{
		// velocità: attenzione alla convenzione dei segni!!
		if ( DOF_v[ 0 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = 1./( 3.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = -1./( 3.0 * s );
		}
		if ( DOF_v[ 1 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = 1./( 3.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = -1./( 3.0 * s );
		}
		if ( DOF_v[ 2 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = 1./( 3.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = -1./( 3.0 * s );
		}

		// pressione
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = -1./3.;
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = -1./3.;
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = -1./3.;

		// pressione media
		(*Aup_i) ( 0 , Index ) = 1.;
	}

	return;

} //setAup_i


void setAup_i ( sparseMatrixPtr_Type& Aup_i,
				size_type id, sizeVector_Type& ID,
				scalarVector_Type& DOF, scalarVector_Type& DOF_v,
				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,
				Matrix4d T, const size_type Index,
				scalar_type s )
{
	size_type id0, id1, id2, id3;

	id0 = ID[0];
	id1 = ID[1];
	id2 = ID[2];
	id3 = ID[3];

	if( id != 4 )
	{
		(*Aup_i) ( 0 , shiftIntersect[ id ] + DOF_v[ id ] )  = 1.;
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( id , 0 );
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( id , 1 );
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = 1.*T( id , 2 );
		(*Aup_i) ( 0 , shiftIntersect[ id3 ] + DOF[ 3 ] + fractureNumberGlobalDOFVelocity [ id3 ] ) = 1.*T( id , 3 );
		(*Aup_i) ( 0 , Index ) = -1.*( T( id , 0 ) + T( id , 1 ) + T( id , 2 ) + T( id, 3 ));
	}
	else
	{
		// velocità: attenzione alla convenzione dei segni!!
		if ( DOF_v[ 0 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v[ 0 ] )  = -1./( 4.0 * s );
		}
		if ( DOF_v[ 1 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v[ 1 ] )  = -1./( 4.0 * s );
		}
		if ( DOF_v[ 2 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF_v[ 2 ] )  = -1./( 4.0 * s );
		}
		if ( DOF_v[ 3 ] == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id3 ] + DOF_v[ 3 ] )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id3 ] + DOF_v[ 3 ] )  = -1./( 4.0 * s );
		}

		// pressione
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = -1./4.;
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF[ 1 ] + fractureNumberGlobalDOFVelocity [ id1 ] ) = -1./4.;
		(*Aup_i) ( 0 , shiftIntersect[ id2 ] + DOF[ 2 ] + fractureNumberGlobalDOFVelocity [ id2 ] ) = -1./4.;
		(*Aup_i) ( 0 , shiftIntersect[ id3 ] + DOF[ 3 ] + fractureNumberGlobalDOFVelocity [ id3 ] ) = -1./4.;

		// pressione media
		(*Aup_i) ( 0 , Index ) = 1.;

	}

	return;

} //setAup_i


void setAup_i ( sparseMatrixPtr_Type& Aup_i, size_type id,
				FracturePtrContainer_Type& Fracture,
				scalarVector_Type& DOF_p0, scalarVector_Type& DOF_v0,
				scalar_type& DOF_p1, scalar_type& DOF_v1,
				sizeVector_Type& shiftIntersect, sizeVector_Type& fractureNumberGlobalDOFVelocity,
				const Matrix4d& T, const size_type Index,
				scalar_type s )
{

	size_type id0 = Fracture[ 0 ]->getID();
	size_type id1 = Fracture[ 1 ]->getID();

	if( id != 4 )
	{
		if( id == 0  || id == 2 )
		{
			if( id == 0)
			{
				(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 0 ] )  = 0.5;
				(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 3 ] )  = 0.5;
			}
			else
			{
				(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 1 ] )  = 0.5;
				(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 2 ] )  = 0.5;
			}
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( id , 0 );
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( id , 1 );
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( id , 2 );
			(*Aup_i) ( 0 , Index + 1 ) = 1.*T( id , 3 );
			(*Aup_i) ( 0 , Index ) = -1.*( T( id , 0 ) + T( id , 1 ) + T( id , 2 ) + T( id , 3 ) );
		}
		else
		{
			if ( id == 1 )
			{
				(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = 1.;
				//nel caso id == 3 imponiamo condizioni no-slip alla parete
			}
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = 1.*T( id , 0 );
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = 1.*T( id , 1 );
			(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = 1.*T( id , 2 );
			(*Aup_i) ( 0 , Index + 1 ) = 1.*T( id , 3 );
			(*Aup_i) ( 0 , Index ) = -1.*( T( id , 0 ) + T( id , 1 ) + T( id , 2 ) + T( id , 3 ) );
		}

	}
	else
	{
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 0 ] )  = -1./( 8.0 * s );
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 3 ] )  = 1./( 8.0 * s );
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 1 ]  )  = 1./( 8.0 * s );
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_v0[ 2 ]  )  = -1./( 8.0 * s );

		if ( DOF_v1 == 0 )
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = 1./( 4.0 * s );
		}
		else
		{
			(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_v1 )  = -1./( 4.0 * s );
		}


		// pressione
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 0 ] + fractureNumberGlobalDOFVelocity [ id0 ] ) = -1./4.;
		(*Aup_i) ( 0 , shiftIntersect[ id1 ] + DOF_p1 + fractureNumberGlobalDOFVelocity [ id1 ] ) = -1./4.;
		(*Aup_i) ( 0 , shiftIntersect[ id0 ] + DOF_p0[ 1 ] ) = -1./4.;
		(*Aup_i) ( 0 , Index + 1 ) = -1./4.;

		// pressione media
		(*Aup_i) ( 0 , Index ) = 1.;
	}

	return;

} //setAup_i


void darcy_dataF ( scalarVectorPtr_Type &Bstress,
                   scalarVectorPtr_Type &Bvel,
                   const BCHandlerPtr_Type& bcHandler,
                   const FractureHandlerPtr_Type& fracture,
                   const scalar_type& gammaU,
                   const scalar_type& invK,
                   const scalarVectorPtr_Type& Pneumann,
                   const scalarVectorPtr_Type& v_diri )
{

    // ----------------- Penalty    ---------------

    const scalar_type fractureID = fracture->getID();
    const size_type shiftVelocity = fracture->getMeshFEM2().nb_dof();
    const size_type shiftCoefficinents = fracture->getMeshFEM().nb_dof();

    scalarVector_Type etaGammaUinvh(shiftCoefficinents, 0.), Bvel_tot( shiftVelocity, 0.);

    for ( size_type i = 0; i < shiftCoefficinents; ++i )
    {
        etaGammaUinvh [ i ] = invK * gammaU * fracture->getInverseMeshSize(i);
    }

    getfem::generic_assembly assem2;

    assem2.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    // Assign the M_mediumMesh integration method
    assem2.push_mi(fracture->getIntegrationMethod2());

    // Assign the M_mediumMesh finite element space
    assem2.push_mf(fracture->getMeshFEM2());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem2.push_mf(fracture->getMeshFEM());
    assem2.push_mf(fracture->getMeshFEM());

    // Assign the coefficients
    assem2.push_data(etaGammaUinvh);
    assem2.push_data(*v_diri);

    // Set the matrices to save the evaluations
    assem2.push_vec(Bvel_tot);

    // Assemble in each sub region
    const size_type shiftDirichlet = bcHandler->getFractureBC(fractureID)->getDirichlet().size();
    for ( size_type bndID = 0; bndID < shiftDirichlet; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getDirichlet(bndID);
        assem2.assembly( fracture->getMeshFEM2().linked_mesh().get_mpi_sub_region( val));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bvel) [ i ] += Bvel_tot [ i ];
    }

    cout << "DARCY :: DATA (penal. bound.)    [OK]" << endl;

    // ----------------- Ext Stress ---------------

    const size_type shiftMapFactor1 = fracture->getMagnificationMapFactor1().size();
    const size_type shiftMapFactor2 = fracture->getMagnificationMapFactor2().size();

    scalarVector_Type Bs( shiftVelocity, 0. ), coefx( shiftMapFactor1 );
    scalarVector_Type coefy( shiftMapFactor2 );

    getfem::generic_assembly assemb;

    if ( fracture->getMeshFEM2().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Base(#2))"
            "(:,1, h).p(h)");
    }
    else
    {

        //attenzione, quando integro sulla frattura devo tenere conto che quella vera non è flat quindi le lunghezze/aree devono essere convertite
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
            coefy [ i ] = 1 / fracture->getMagnificationMapFactor2(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
            "(:,k, k, h).p(h)");
    }

    // Assign the M_mediumMesh integration method
    assemb.push_mi(fracture->getIntegrationMethod2());

    // Assign the M_mediumMesh finite element space
    assemb.push_mf(fracture->getMeshFEM2());

    // Assign the M_mediumMesh finite element space for the coefficients
    assemb.push_mf(bcHandler->getFractureBC(fracture->getID())->getMeshFEM());

    // Assign the coefficients
    assemb.push_data(*Pneumann);

    // Set the vector to save the evaluations
    assemb.push_vec(Bs);

    // Assemble in each sub region
    const size_type shiftNeumann = bcHandler->getFractureBC(fractureID)->getNeumann().size();

    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getNeumann( bndID);
        assemb.assembly( fracture->getMeshFEM2().linked_mesh().get_mpi_sub_region( val ));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bstress) [ i ] = Bs [ i ];
    }

    return;

} // darcy_dataF


//termine sorgente per la frattura
void assembling_Source_BoundaryF ( scalarVectorPtr_Type& D,
                                   const scalarVectorPtr_Type& source,
                                   const FractureHandlerPtr_Type& fracture,
                                   const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = fracture->getMeshFEM().nb_dof();
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();

    scalarVector_Type D_(shiftPressure, 0.0), invF(shiftMapFactor, 0);

    generic_assembly assem_Source, assem_Vx, assem_Vy;

    if ( fracture->getMeshFEM().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i));
        }

        /*
         * assembla il termine noto
         *
         * 		sum(i) { - (f_i, q_i) }
         */
        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i) * fracture->getMagnificationMapFactor2(i));
        }

        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }

    // Assign the M_mediumMesh integration method
    assem_Source.push_mi(fracture->getIntegrationMethod());

    // Assign the M_mediumMesh finite element space
    assem_Source.push_mf(fracture->getMeshFEM());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_Source.push_mf(fracture->getMeshFEM());
    assem_Source.push_mf(fracture->getMeshFEM());

    // Assign the coefficients
    assem_Source.push_data(*source);
    assem_Source.push_data(invF);

    // Set the vector to save the evaluations
    assem_Source.push_vec(D_);

    // Computes the matrices
    assem_Source.assembly ( uncutRegionFlag );

    gmm::add(D_, gmm::sub_vector(*D, gmm::sub_interval(0, shiftPressure)));

    return;

} // assembling_Source_BoundaryF


}// namespace getfem
