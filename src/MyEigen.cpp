// MyEigen.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "MyEigen.h"

//**************************************************************************
// MyEigen class.

void MyEigen::flowVelocity4D(Lattice *lat, Group *group, Parameters *param, int it)
{
  int N = param->getSize();
  int pos;
  int changeSign;
  int count;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double x, y, ux=0., uy=0., ueta=0., utau=0.;
  double dtau = param->getdtau();
  gsl_complex square;
  gsl_complex factor;
  gsl_complex euklidiansquare;
  double eps;
  gsl_complex z_aux;
  gsl_complex tau2;
  int foundU;
  // stringstream stru_name;
  // stru_name << "u" << param->getMPIRank() << ".dat";
  // string u_name;
  // u_name = stru_name.str();
  // ofstream fout(u_name.c_str(),ios::app);
  // stringstream strux_name;
  // strux_name << "u-x" << param->getMPIRank() << ".dat";
  // string ux_name;
  // ux_name = strux_name.str();
  // ofstream foutx(ux_name.c_str(),ios::app); 
  double averageux=0.;
  double averageuy=0.;
  double averageueta=0.;
  double averageeps=0.;
  count = 0;
  double uxd, uyd, uetad;
  
  for (int si=0; si<N; si++)
    {
      for (int sj=0; sj<N; sj++)
	{

	  x = -L/2.+a*si;
	  y = -L/2.+a*sj;
	
	  //	  if(si%10==0 && sj%10==0)
	  //  cout << si << " " << sj << endl;

	  pos = si*N+sj;
	  GSL_SET_COMPLEX(&square, 0, 0);





	  double data[] = { lat->cells[pos]->getTtautau(), -lat->cells[pos]->getTtaux(), -lat->cells[pos]->getTtauy(), -(it*dtau*a)*(it*dtau*a)*lat->cells[pos]->getTtaueta(),
	  		    lat->cells[pos]->getTtaux()  , -lat->cells[pos]->getTxx()  , -lat->cells[pos]->getTxy()  , -(it*dtau*a)*(it*dtau*a)*lat->cells[pos]->getTxeta(), 
	  		    lat->cells[pos]->getTtauy()  , -lat->cells[pos]->getTxy()  , -lat->cells[pos]->getTyy()  , -(it*dtau*a)*(it*dtau*a)*lat->cells[pos]->getTyeta(),
	  		    lat->cells[pos]->getTtaueta(), -lat->cells[pos]->getTxeta(), -lat->cells[pos]->getTyeta(), -(it*dtau*a)*(it*dtau*a)*lat->cells[pos]->getTetaeta()
	  };


	  
	  gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4); //matrix
	  
	  gsl_vector_complex *eval = gsl_vector_complex_alloc (4); //eigenvalues are components of this vector
	  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4); //eigenvectors are columns of this matrix
	  
	  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4); //workspace
	  
	  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w); // solve for eigenvalues and eigenvectors (without 'v' only compute eigenvalues)
	  
	  gsl_eigen_nonsymmv_free (w); // free memory associated with workspace
	  

	  //set to 'zero' 
	  lat->cells[pos]->setEpsilon( lat->cells[pos]->getTtautau() );
	  lat->cells[pos]->setutau(1.);
	  lat->cells[pos]->setux(0.);
	  lat->cells[pos]->setuy(0.);
	  lat->cells[pos]->setueta(0.);

	  //gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);  // sort by ascending order in magnitude
	  
	  { // output:
	    int i;
	    foundU=0;
	    for (i = 0; i < 4; i++)
	      {
		gsl_complex eval_i = gsl_vector_complex_get (eval, i);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
		
		GSL_SET_COMPLEX(&square, 0, 0);
		GSL_SET_COMPLEX(&tau2, a*it*dtau*a*it*dtau, 0);
		GSL_SET_COMPLEX(&euklidiansquare, 0, 0);
		
		for (int j = 0; j < 4; ++j)
		  {	
		    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
		    z_aux = gsl_complex_mul(tau2,z);	
		    euklidiansquare = gsl_complex_add(euklidiansquare, gsl_complex_mul(z,z));
		    
		    if (j==0)
		      square = gsl_complex_add(square, gsl_complex_mul(z,z));
		    else if (j<3)
		      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
		    else
		      square = gsl_complex_sub(square, gsl_complex_mul(z_aux,z));
		  }
		
		GSL_SET_COMPLEX(&factor, sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square))), 0);
		if(GSL_REAL(square)>0 )
		  {
		    eps = GSL_REAL(eval_i);///sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square)));
		    if(abs(GSL_IMAG(eval_i))>0.001 && si>0 && sj>0 && si<N-5 && sj<N-5)
		      {
			// cout << " ********* Warning: imaginary energy density" << endl;
			// cout << "Re=" << GSL_REAL(eval_i) << ", Im=" << GSL_IMAG(eval_i) << endl;
			// cout << si << " " << sj << endl << endl;
			// cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() 
			//      << " " << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTtaueta() << endl;
			// cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() 
			//      << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTxeta() << endl;
			// cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() 
			//      << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta() << endl;
			// cout << lat->cells[pos]->getTtaueta() << " " << lat->cells[pos]->getTxeta() 
			//      << " " << lat->cells[pos]->getTyeta() << " " << lat->cells[pos]->getTetaeta() << endl;
			eps = lat->cells[pos]->getTtautau();
		      }
		 
		 //    if(si==180 && sj==109)
// 		      {
// 			cout << si << " " << sj << endl;
// 			cout << "Re=" << GSL_REAL(eval_i) << ", Im=" << GSL_IMAG(eval_i) << endl;
// 			cout << si << " " << sj << endl << endl;
// 			cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() << endl;
// 			cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() << endl;
// 			cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTyy() << endl;
 	   
// 		      }
		  }
				
		GSL_SET_COMPLEX(&square, 0, 0);
	
		for (int j = 0; j < 4; ++j)
		  {	
		    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
		    z = gsl_complex_mul(z,factor);
		    z_aux = gsl_complex_mul(tau2,z);	
		   
		    if (j==0)
		      square = gsl_complex_add(square, gsl_complex_mul(z,z));
		    else if (j<3)
		      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
		    else
		      square = gsl_complex_sub(square, gsl_complex_mul(z_aux,z));
		  }
		
	// 	cout<< "i=" << i << ", square=" << GSL_REAL(square) << ", vector=" << endl;

// 		for (int j = 0; j < 4; ++j)
// 		  {	
// 		    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
// 		    cout << GSL_REAL(z) << " " << GSL_IMAG(z) << endl;
// 		  }

		changeSign=0;
		// for the time-like eigenvector do the following (this is the flow velocity)
		if(GSL_REAL(square)>0)
		  {
		    foundU+=1;
		    for (int j = 0; j < 4; ++j)
		      {	
			gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
			z = gsl_complex_mul(z,factor);
			
			if(j==0 && GSL_REAL(z)<0)
			  {
			    changeSign=1;
			    GSL_SET_COMPLEX(&z, -1.*GSL_REAL(z), -1.*GSL_IMAG(z));
			  }
			
			if(j>0 && changeSign==1)
			  {
			    GSL_SET_COMPLEX(&z, -1.*GSL_REAL(z), -1.*GSL_IMAG(z));
			  }
			
			if(j==0)
			  {
			    if(eps>0.1)
			      utau = GSL_REAL(z);
			    else 
			      utau=1.;
			  }
			if(j==1)
			  {
			    if(eps>0.1)
			      ux = GSL_REAL(z);
			    else 
			      ux=0.;
			  }
			if(j==2)
			  {
			    if(eps>0.1)
			      uy = GSL_REAL(z);
			    else 
			      uy=0.;
			  }
			if(j==3)
			  {
			    if(eps>0.1)
			      ueta = GSL_REAL(z);
			    else 
			      ueta=0.;
			  }
		
	
			//	if(si%10==0 && sj%10==0)
			// printf("%g\n", GSL_REAL(z));
			
			if(abs(GSL_IMAG(z))>0.001 && eps>0.001 && si>10 && sj>10 && si<N-10  && sj<N-10)
			  {
			    // cout << si << " " << sj << endl << endl;
			    // cout << "u(" << j << ")=" << GSL_REAL(z) << ", " << GSL_IMAG(z) << endl;
			    // printf("u_mu u^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));
			    utau=1.;
			    ux = 0.;
			    uy = 0.;
			    ueta=0.;
	      
			    eps = lat->cells[pos]->getTtautau();
			  }
		
	
	// 		if (GSL_IMAG(z)!=0)
// 			  {
// 			    cout << "warning - imaginary flow" << endl;
// 			    cout << GSL_IMAG(z) << endl;
// 			  }
		      }
		    
		    //cout << it*dtau << endl;
		    //    cout << "it=" << it << endl;
	
		    //		      cout <<  x << " " << y << " " << eps << " " << utau << " " << ux << " " << uy << " " << a*it*dtau*ueta << endl;

		    // double tau = a*it*dtau;
		    // uxd = -ux;
		    // uyd = -uy;
		    // uetad = -tau*tau*ueta;

		    // double eps2 = utau * utau * lat->cells[pos]->getTtautau() + 2.* utau * uxd * lat->cells[pos]->getTtaux()
		    //   + 2.* utau * uyd * lat->cells[pos]->getTtauy() + 2.* utau *uetad * lat->cells[pos]->getTtaueta()/a
		    //   + uxd * uxd * lat->cells[pos]->getTxx() + 2. * uxd * uyd * lat->cells[pos]->getTxy() 
		    //   + 2. * uxd * uetad * lat->cells[pos]->getTxeta()/a + uyd * uyd * lat->cells[pos]->getTyy() 
		    //   + 2. * uyd * uetad * lat->cells[pos]->getTyeta()/a + uetad *uetad * lat->cells[pos]->getTetaeta()/a/a;
		    
		    // cout << eps << " " << eps2 << " " << utau << " " << ux << " " << uy << " " << a*it*dtau*ueta << endl;

		    lat->cells[pos]->setutau(utau);
		    lat->cells[pos]->setux(ux);
		    lat->cells[pos]->setuy(uy);
		    lat->cells[pos]->setueta(ueta);
		    lat->cells[pos]->setEpsilon(eps);

		    averageux+=ux*ux*eps;
		    averageuy+=uy*uy*eps;
		    averageueta+=ueta*ueta*eps*it*dtau*a*it*dtau*a;
		    averageeps+=eps;

		    count ++;
		    // if(si>0 && sj>0 && param->getWriteOutputs()>0 && foundU==1) // 
		    //   {
		    // 	fout << x << " " << y << " " << eps << " " << utau << " " << ux << " " << uy << " " << a*it*dtau*ueta << " " 
		    // 	     << " " << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() 
		    // 	     << " " << lat->cells[pos]->getTtaueta()/a << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() 
		    // 	     << " " << lat->cells[pos]->getTxeta()/a << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta()/a << " " << lat->cells[pos]->getTetaeta()/a/a << endl;
		      }
		    // if(sj == N/2 && param->getWriteOutputs()>0)
		    //   {
		    // 	foutx << a*it*dtau << " " << x << " " << y << " " << ux << " " << uy << " " << a*it*dtau*ueta << " " 
		    // 	      << eps << " " << lat->cells[pos]->getTtautau() << " " << " " << lat->cells[pos]->getTetaeta() << endl;
		    // 	if(si == N-1)
		    // 	  {
		    // 	    foutx << endl;
		    // 	  }			    
		    //   }
		    
	      }

	  }

	    // write Tmunu in case no u was found
	    if ( foundU==0 ) 
	      {
	      // 	if(si>0 && sj>0 && param->getWriteOutputs()>0) // 
	      // 	  {
	      // 	    fout << x << " " << y << " " << eps << " " << utau << " " << ux << " " << uy << " " << a*it*dtau*ueta << " " 
	      // 		 << " " << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() 
	      // 		 << " " << lat->cells[pos]->getTtaueta()/a << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() 
	      // 		 << " " << lat->cells[pos]->getTxeta()/a << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta()/a << " " << lat->cells[pos]->getTetaeta()/a/a << endl;
	      // 	  }
	      // }

	    if(si==N/2 && sj==N/2)
	      {	cout << si << " " << sj << endl << endl;
		cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() 
		     << " " << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTtaueta() << endl;
		cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() 
		     << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTxeta() << endl;
		cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() 
		     << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta() << endl;
		cout << lat->cells[pos]->getTtaueta() << " " << lat->cells[pos]->getTxeta() 
		     << " " << lat->cells[pos]->getTyeta() << " " << lat->cells[pos]->getTetaeta() << endl;
		cout << "ux=" << ux << endl;
		cout << "uy=" << uy << endl;
		cout << "ueta=" << ueta << endl;
	      }
	    
	    // gsl_vector_complex_free (eval);
	    // gsl_matrix_complex_free (evec);
	    
	      }
	  
	  // if(foundU==0)
	  //   {
	  //     if(si>10 && sj>10 && si<N-10  && sj<N-10)
	  // 	cout << "Did not find a flow velocity (no time-like eigenvector)... in cell " << si << ", " << sj << endl;
	  //     if(si%6==0 && sj%6==0 && si>0 && sj>0)
	  // 	{
	  // 	  fout << a*it*dtau << " " << x << " " << y << " " << 0. << " " << 0. << " " << 0. << " " 
	  // 	       << 0. << " " << 0. << " " << 0. << endl;
	  // 	}
	      
	  //   }

	  // if(foundU>1 && si>10 && sj>10 && si<N-10 && sj<N-10)
	  //   {
	  //     //	      cout << "Found more than one flow velocity (more than one time-like eigenvector)... in cell " << si << ", " << sj << endl;
	      
	  //   }
	  
	  
	    // compute pi^{\mu\nu}
	    if(utau==1 && ux ==0 && uy == 0 && ueta == 0)
	      {
		lat->cells[pos]->setpitautau(0.);
		lat->cells[pos]->setpixx(0.);
		lat->cells[pos]->setpiyy(0.);
		lat->cells[pos]->setpietaeta(0.);
		
		lat->cells[pos]->setpitaux(0.);
		lat->cells[pos]->setpitauy(0.);
		lat->cells[pos]->setpitaueta(0.);
		
		lat->cells[pos]->setpixeta(0.);
		lat->cells[pos]->setpixy(0.);
		lat->cells[pos]->setpiyeta(0.);
	      }
	    else
	      {
		lat->cells[pos]->setpitautau(lat->cells[pos]->getTtautau() - 4./3.*eps*utau*utau + eps/3.);
		lat->cells[pos]->setpixx(lat->cells[pos]->getTxx() - 4./3.*eps*ux*ux - eps/3.);
		lat->cells[pos]->setpiyy(lat->cells[pos]->getTyy() - 4./3.*eps*uy*uy - eps/3.);
		lat->cells[pos]->setpietaeta(lat->cells[pos]->getTetaeta() - 4./3.*eps*ueta*ueta - eps/3./it/dtau/a/it/dtau/a);
		
		lat->cells[pos]->setpitaux(lat->cells[pos]->getTtaux() - 4./3.*eps*utau*ux);
		lat->cells[pos]->setpitauy(lat->cells[pos]->getTtauy() - 4./3.*eps*utau*uy);
		lat->cells[pos]->setpitaueta(lat->cells[pos]->getTtaueta() - 4./3.*eps*utau*ueta);
		
		lat->cells[pos]->setpixeta(lat->cells[pos]->getTxeta() - 4./3.*eps*ux*ueta);
		lat->cells[pos]->setpixy(lat->cells[pos]->getTxy() - 4./3.*eps*ux*uy);
		lat->cells[pos]->setpiyeta(lat->cells[pos]->getTyeta() - 4./3.*eps*uy*ueta);
	      }
	    

	    

	    // cout << "pi diag: " << lat->cells[pos]->getpitautau() << " " << lat->cells[pos]->getpixx() << " " 
	    // 	 << lat->cells[pos]->getpiyy() << " " << lat->cells[pos]->getpietaeta() << endl;
	    // cout << "T diag: " << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTxx() << " " 
	    // 	 << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTetaeta() << endl;
	    // cout << "pi trace: " << lat->cells[pos]->getpitautau() - lat->cells[pos]->getpixx() - lat->cells[pos]->getpiyy() 
	    //   -it*dtau*a*it*dtau*a*lat->cells[pos]->getpietaeta() << endl;
	    
	 
	    gsl_vector_complex_free (eval);
	    gsl_matrix_complex_free (evec);
	
	}
    }

  cout << it*dtau*a << " average u^x=" << sqrt(averageux/averageeps) << endl;
  cout << it*dtau*a << " average u^y=" << sqrt(averageuy/averageeps) << endl;
  cout << it*dtau*a << " average tau u^eta=" << sqrt(averageueta/averageeps) << endl;

  //  fout.close();
  //  foutx.close();

  double maxtime = param->getMaxtime(); // maxtime is in fm
  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));

  // output for hydro
  //  if(it==itmax && param->getWriteOutputs() > 0)
  if(param->getWriteOutputs() > 0)
    {
      double alphas = param->getalphas();
      double g = param->getg();
      double gfactor;
      double Qs;
      int hx = param->getSizeOutput();
      int hy = hx;
      int heta = param->getEtaSizeOutput();
      double hL = param->getLOutput();
      double deta = param->getDetaOutput();
      double c = param->getc();
      double muZero = param->getMuZero();
      double g2mu2A, g2mu2B;
      double PI = param->getPi();

      if(hL>L)
	cout << "WARNING: hydro grid length larger than the computed one." << endl;
      
      int xpos, ypos, xposUp, yposUp, pos1, pos2, pos3;
      double fracx, fracy, x1, x2;
      double xlow, xhigh, ylow, yhigh;
      int pos4;
      double resultE, resultutau, resultux, resultuy, resultueta;
      double resultpi00, resultpi0x, resultpi0y, resultpi0eta;
      double resultpixy, resultpixeta, resultpiyeta, resultpixx, resultpiyy, resultpietaeta;

      double ha;
      ha = hL/static_cast<double>(hx);
      
      stringstream streuH_name;
      streuH_name << "epsilon-u-Hydro-t" << it*dtau*a << "-" << param->getMPIRank() << ".dat";
      string euH_name;
      euH_name = streuH_name.str();

      ofstream foutEps2(euH_name.c_str(),ios::out); 
      foutEps2 << "# dummy " << 1 << " etamax= " << heta
	       << " xmax= " << hx << " ymax= " << hy << " deta= " << deta 
	       << " dx= " << ha << " dy= " << ha << endl; 
      
       for(int ieta=0; ieta<heta; ieta++) // loop over all positions
	{
	  for(int ix=0; ix<hx; ix++) // loop over all positions
	    {
	      for(int iy=0; iy<hy; iy++)
		{
		  x = -hL/2.+ha*ix;
		  y = -hL/2.+ha*iy;
		  
		  //		  cout << ix << " " << iy << endl;
		  
		  if (abs(x) < L/2. && abs(y) < L/2.)
		    {
		      xpos = static_cast<int>(floor((x+L/2.)/a+0.0000000001));
		      ypos = static_cast<int>(floor((y+L/2.)/a+0.0000000001));
		      
		      
		      if(xpos<N-1)
			xposUp = xpos+1;
		      else
			xposUp = xpos;
		      
		      if(ypos<N-1)
			yposUp = ypos+1;
		      else
			yposUp = ypos;
		      
		      xlow = -L/2.+a*xpos;
		      ylow = -L/2.+a*ypos;
		      
		      xhigh = -L/2.+a*xposUp;
		      yhigh = -L/2.+a*yposUp;
		      
		      fracx = (x-xlow)/ha;
		      
		      pos1 = xpos*N+ypos;
		      pos2 = xposUp*N+ypos;
		      pos3 = xpos*N+yposUp;
		      pos4 = xposUp*N+yposUp;
		      

		      // -----------------------------epsilon--------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getEpsilon());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getEpsilon());
		      else
			x2 = 0.;
		      
		      fracy = (y-ylow)/ha;
		      
		      resultE = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------utau------------------------------------ //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getutau())+fracx*(lat->cells[pos2]->getutau());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getutau())+fracx*(lat->cells[pos4]->getutau());
		      else
			x2 = 0.;
		      
		      resultutau = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------ux-------------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getux())+fracx*(lat->cells[pos2]->getux());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getux())+fracx*(lat->cells[pos4]->getux());
		      else
			x2 = 0.;
		      
		      resultux = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------uy-------------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getuy())+fracx*(lat->cells[pos2]->getuy());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getuy())+fracx*(lat->cells[pos4]->getuy());
		      else
			x2 = 0.;
		      
		      resultuy = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------ueta------------------------------------ //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getueta())+fracx*(lat->cells[pos2]->getueta());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getueta())+fracx*(lat->cells[pos4]->getueta());
		      else
			x2 = 0.;
		      
		      resultueta = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pitautau-------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpitautau())+fracx*(lat->cells[pos2]->getpitautau());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpitautau())+fracx*(lat->cells[pos4]->getpitautau());
		      else
			x2 = 0.;
		      
		      resultpi00 = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pitaux---------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpitaux())+fracx*(lat->cells[pos2]->getpitaux());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpitaux())+fracx*(lat->cells[pos4]->getpitaux());
		      else
			x2 = 0.;
		      
		      resultpi0x = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pitauy---------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpitauy())+fracx*(lat->cells[pos2]->getpitauy());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpitauy())+fracx*(lat->cells[pos4]->getpitauy());
		      else
			x2 = 0.;
		      
		      resultpi0y = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pitaueta-------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpitaueta())+fracx*(lat->cells[pos2]->getpitaueta());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpitaueta())+fracx*(lat->cells[pos4]->getpitaueta());
		      else
			x2 = 0.;
		      
		      resultpi0eta = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pixy---------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpixy())+fracx*(lat->cells[pos2]->getpixy());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpixy())+fracx*(lat->cells[pos4]->getpixy());
		      else
			x2 = 0.;
		      
		      resultpixy = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pixeta---------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpixeta())+fracx*(lat->cells[pos2]->getpixeta());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpixeta())+fracx*(lat->cells[pos4]->getpixeta());
		      else
			x2 = 0.;
		      
		      resultpixeta = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------piyeta---------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpiyeta())+fracx*(lat->cells[pos2]->getpiyeta());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpiyeta())+fracx*(lat->cells[pos4]->getpiyeta());
		      else
			x2 = 0.;
		      
		      resultpiyeta = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pixx------------------------------------ //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpixx())+fracx*(lat->cells[pos2]->getpixx());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpixx())+fracx*(lat->cells[pos4]->getpixx());
		      else
			x2 = 0.;
		      
		      resultpixx = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------piyy------------------------------------ //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpiyy())+fracx*(lat->cells[pos2]->getpiyy());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpiyy())+fracx*(lat->cells[pos4]->getpiyy());
		      else
			x2 = 0.;
		      
		      resultpiyy = (1.-fracy)*x1+fracy*x2;
		      
		      // -----------------------------pietaeta-------------------------------- //
		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
			x1 = (1-fracx)*(lat->cells[pos1]->getpietaeta())+fracx*(lat->cells[pos2]->getpietaeta());
		      else
			x1 = 0.;
		      
		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
			x2 = (1-fracx)*(lat->cells[pos3]->getpietaeta())+fracx*(lat->cells[pos4]->getpietaeta());
		      else
			x2 = 0.;
		      
		      resultpietaeta = (1.-fracy)*x1+fracy*x2;
		      
		      if(resultutau<1.)
			{
			  resultutau=1.;
			  resultux=0.;
			  resultuy=0;
			  resultueta=0;
			}
		      
		      
		      if(param->getRunningCoupling())
			{
			  
			  // if(pos>0 && pos<(N-1)*N+N-1)
			  // 	{
			  // 	  g2mu2A = lat->cells[pos]->getg2mu2A();
			  // 	}
			  // else 
			  // 	g2mu2A = 0;
			  
			  // if(pos>0 && pos<(N-1)*N+N-1)
			  // 	{
			  // 	  g2mu2B = lat->cells[pos]->getg2mu2B();
			  // 	}
			  // else
			  // 	g2mu2B = 0;
			  
			  // if(param->getRunWithQs()==2)
			  // 	{
			  // 	  if(g2mu2A > g2mu2B)
			  // 	    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  // 	  else
			  // 	    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  // 	}
			  // else if(param->getRunWithQs()==0)
			  // 	{
			  // 	  if(g2mu2A < g2mu2B)
			  // 	    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  // 	  else
			  // 	    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  // 	}
			  // else if(param->getRunWithQs()==1)
			  // 	{
			  // 	  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  // 	}
			  
			  // 3 flavors
			  //	      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
			  //gfactor = g*g/(4.*PI*alphas);
			  
			  
			  // run with average Q_s only ! local makes no sense here (stuff has moved in the mean time)
			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));	  
			  gfactor = g*g/(4.*PI*alphas);
			  //     cout << "Qs=" << Qs << endl;
			  //cout << alphas << endl;
			}
		      else
			gfactor = 1.;
		      
		     	      
		      if(abs(0.1973269718*resultE*gfactor) > 0.0000000001)
			foutEps2 << -(heta-1)/2.*deta+deta*ieta << " " << x << " " << y << " " 
				 << abs(0.1973269718*resultE*gfactor) 
				 << " " << resultutau << " " << resultux << " " << resultuy << " " << resultueta 
				 << " " << resultpi00*gfactor << " " << resultpi0x*gfactor
				 << " " << resultpi0y*gfactor << " " << resultpi0eta*gfactor 
				 << " " << resultpixx*gfactor << " " << resultpixy*gfactor 
				 << " " << resultpixeta*gfactor << " " << resultpiyy*gfactor 
				 << " " << resultpiyeta*gfactor << " " << resultpietaeta*gfactor 
				 << endl; 
		      //      cout << "test" << ix << " " << iy << endl;
		      else
			foutEps2 << -(heta-1)/2.*deta+deta*ieta << " " << x << " " << y << " " 
				 << 0. << " " << 1. << " " << 0. << " " << 0. << " " << 0. 
				 << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. 
				 << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
		    } 
		  else
		    foutEps2 << -(heta-1)/2.*deta+deta*ieta << " " << x << " " << y << " " 
			     << 0. << " " << 1. << " " << 0. << " " << 0. << " " << 0. 
			     << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. 
			     << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
		}
	    }
	  foutEps2 << endl;
	}
       foutEps2.close();
    }

  cout << "Wrote outputs" << endl;
  // done output for hydro
}


void MyEigen::test()
{
  gsl_complex square;
  gsl_complex factor;
  gsl_complex euklidiansquare;
  GSL_SET_COMPLEX(&square, 0, 0);

  double data[] = { 0.1,  -0.001, -0.01, 
		    0.001, -0.1,  -0.0001, 
		    0.01, -0.0001,  -0.1 };

  gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3); //matrix

  gsl_vector_complex *eval = gsl_vector_complex_alloc (3); //eigenvalues are components of this vector
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (3, 3); //eigenvectors are columns of this matrix
  
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (3); //workspace
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w); // solve for eigenvalues and eigenvectors (without 'v' only compute eigenvalues)
  
  gsl_eigen_nonsymmv_free (w); // free memory associated with workspace
  
  //gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);  // sort by ascending order in magnitude
  
  { // output:
    int i;
    
    for (i = 0; i < 3; i++)
      {
	gsl_complex eval_i = gsl_vector_complex_get (eval, i);
	gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

	printf ("before eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
	GSL_SET_COMPLEX(&square, 0, 0);
	GSL_SET_COMPLEX(&euklidiansquare, 0, 0);

	for (int j = 0; j < 3; ++j)
	  {	
	    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	    printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
	    euklidiansquare = gsl_complex_add(euklidiansquare, gsl_complex_mul(z,z));
	    
	    if (j==0)
	      square = gsl_complex_add(square, gsl_complex_mul(z,z));
	    else
	      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
	  }

	if(GSL_REAL(square)>0)
	  printf("before u_mu u^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));
	else 
	  printf("before z_mu z^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));

	GSL_SET_COMPLEX(&factor, sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square))), 0);
	printf ("eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
	
	GSL_SET_COMPLEX(&square, 0, 0);

	for (int j = 0; j < 3; ++j)
	  {	
	    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	    z = gsl_complex_mul(z,factor);

	    printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
    
	    if (j==0)
	      square = gsl_complex_add(square, gsl_complex_mul(z,z));
	    else
	      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
	  }

	if(GSL_REAL(square)>0)
	  printf("u_mu u^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));
	else 
	  printf("z_mu z^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));
	  

      }
  }
  
  gsl_vector_complex_free (eval);
  gsl_matrix_complex_free (evec);
  
  exit(1);
}


//3D version, without eta components

void MyEigen::flowVelocity(Lattice *lat, Group *group, Parameters *param, int it)
{
  int N = param->getSize();
  int pos;
  int changeSign;
  int count;
  double L = param->getL();
  double a = L/static_cast<double>(N); // lattice spacing in fm
  double x, y, utau, ux, uy;
  double dtau = param->getdtau();
  gsl_complex square;
  gsl_complex factor;
  gsl_complex euklidiansquare;
  double eps;

  // stringstream stru_name;
  // stru_name << "u3D" << param->getMPIRank() << ".dat";
  // string u_name;
  // u_name = stru_name.str();

  // ofstream fout(u_name.c_str(),ios::app); 

  // stringstream strux_name;
  // strux_name << "u-x3D" << param->getMPIRank() << ".dat";
  // string ux_name;
  // ux_name = strux_name.str();

  // ofstream foutx(ux_name.c_str(),ios::app); 

  double averageuT=0.;
  count = 0;
  for (int si=0; si<N; si++)
    {
      for (int sj=0; sj<N; sj++)
	{

	  x = -L/2.+a*si;
	  y = -L/2.+a*sj;
	
	  //	  if(si%10==0 && sj%10==0)
	  //  cout << si << " " << sj << endl;

	  pos = si*N+sj;
	  GSL_SET_COMPLEX(&square, 0, 0);

	  double data[] = { lat->cells[pos]->getTtautau() , -lat->cells[pos]->getTtaux() , -lat->cells[pos]->getTtauy(), 
			    lat->cells[pos]->getTtaux()   , -lat->cells[pos]->getTxx()   , -lat->cells[pos]->getTxy(), 
			    lat->cells[pos]->getTtauy()   , -lat->cells[pos]->getTxy()   , -lat->cells[pos]->getTyy() 
	  };
	  

	  // cout << lat->cells[pos]->getTtautau() << " " <<  -lat->cells[pos]->getTtaux() << " " <<  -lat->cells[pos]->getTtauy() << endl;
	  // if (abs(lat->cells[pos]->getTtautau()) <  abs(-lat->cells[pos]->getTtaux()) ||
	  //     abs(lat->cells[pos]->getTtautau()) <  abs(-lat->cells[pos]->getTtauy()) )
	  //   cout << "********************************" << endl;
	    
	  gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3); //matrix
	  
	  gsl_vector_complex *eval = gsl_vector_complex_alloc (3); //eigenvalues are components of this vector
	  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (3, 3); //eigenvectors are columns of this matrix
	  
	  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (3); //workspace
	  
	  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w); // solve for eigenvalues and eigenvectors (without 'v' only compute eigenvalues)
	  
	  gsl_eigen_nonsymmv_free (w); // free memory associated with workspace
	  
	  //gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);  // sort by ascending order in magnitude
	  
	  //set to 'zero' 
	  lat->cells[pos]->setEpsilon( lat->cells[pos]->getTtautau() );
	  lat->cells[pos]->setutau(1.);
	  lat->cells[pos]->setux(0.);
	  lat->cells[pos]->setuy(0.);
	  
		
	  //{ // output:
	    int i;
	    
	    for (i = 0; i < 3; i++)
	      {
		gsl_complex eval_i = gsl_vector_complex_get (eval, i);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
		
		GSL_SET_COMPLEX(&square, 0, 0);
		GSL_SET_COMPLEX(&euklidiansquare, 0, 0);
		
		for (int j = 0; j < 3; ++j)
		  {	
		    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
		    euklidiansquare = gsl_complex_add(euklidiansquare, gsl_complex_mul(z,z));
		    
		    if (j==0)
		      square = gsl_complex_add(square, gsl_complex_mul(z,z));
		    else
		      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
		  }
		
		GSL_SET_COMPLEX(&factor, sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square))), 0);
		if(GSL_REAL(square)>0)
		  {
		    // cout << si << " " << sj << "sqrt(factor)=" << sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square))) << endl;
		    //    cout << "eps(before)=" << eps << endl;
		    eps = GSL_REAL(eval_i);///sqrt(abs(GSL_REAL(euklidiansquare)/GSL_REAL(square)));
		    //cout << "eps(after)=" << eps << endl;
		    if(abs(GSL_IMAG(eval_i))>0.001 && si>0 && sj>0 && si<N-5 && sj<N-5)
		      {
			// cout << " ********* Warning: imaginary energy density" << endl;
			// cout << "Re=" << GSL_REAL(eval_i) << ", Im=" << GSL_IMAG(eval_i) << endl;
			// cout << si << " " << sj << endl << endl;
			// cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() << endl;
			// cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() << endl;
			// cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTyy() << endl;
			eps = lat->cells[pos]->getTtautau();
		      }
		 
		  //   if(si==180 && sj==109)
// 		      {
// 			cout << si << " " << sj << endl;
// 			cout << "Re=" << GSL_REAL(eval_i) << ", Im=" << GSL_IMAG(eval_i) << endl;
// 			cout << si << " " << sj << endl << endl;
// 			cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() << endl;
// 			cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() << endl;
// 			cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTyy() << endl;
 	   
// 		      }
		  }
				
		GSL_SET_COMPLEX(&square, 0, 0);
		
		for (int j = 0; j < 3; ++j)
		  {	
		    gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
		    z = gsl_complex_mul(z,factor);
		    
		    if (j==0)
		      square = gsl_complex_add(square, gsl_complex_mul(z,z));
		    else
		      square = gsl_complex_sub(square, gsl_complex_mul(z,z));
		  }
		
		//		cout << si << " " << sj << " " << i << " " << GSL_REAL(square) << " " << GSL_REAL(factor) << endl;
	
		changeSign=0;
		// for the time-like eigenvector do the following (this is the flow velocity)
		if(GSL_REAL(square)>0)
		  {
		    for (int j = 0; j < 3; ++j)
		      {	
			gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
			z = gsl_complex_mul(z,factor);
			
			if(j==0 && GSL_REAL(z)<0)
			  {
			    changeSign=1;
			    GSL_SET_COMPLEX(&z, -1.*GSL_REAL(z), -1.*GSL_IMAG(z));
			  }
			
			if(j>0 && changeSign==1)
			  {
			    GSL_SET_COMPLEX(&z, -1.*GSL_REAL(z), -1.*GSL_IMAG(z));
			  }
			
			if(j==0)
			  {
			    if(eps>0.0001)
			      utau = GSL_REAL(z);
			    else 
			      utau=1.;
			  }
			if(j==1)
			  {
			    if(eps>0.0001)
			      ux = GSL_REAL(z);
			    else 
			      ux=0.;
			  }
			if(j==2)
			  {
			    if(eps>0.0001)
			      uy = GSL_REAL(z);
			    else 
			      uy=0.;
			  }
		
	
			//	if(si%10==0 && sj%10==0)
			// printf("%g\n", GSL_REAL(z));
			
			if(abs(GSL_IMAG(z))>0.001 && eps>0.001 && si>0 && sj>0 && si<N-5 && sj<N-5)
			  {
			    // cout << si << " " << sj << endl << endl;
			    // cout << "u(" << j << ")=" << GSL_REAL(z) << ", " << GSL_IMAG(z) << endl;
			    // printf("u_mu u^mu = %g + %gi\n", GSL_REAL(square), GSL_IMAG(square));
			    utau = 1.;
			    ux = 0.;
			    uy = 0.;
			    eps = lat->cells[pos]->getTtautau();
			  }
		
			if(abs(GSL_IMAG(z))>0.001)
			  {
			    utau = 1.;
			    ux = 0.;
			    uy = 0.;
			    eps = lat->cells[pos]->getTtautau();
			  }
		
	// 		if (GSL_IMAG(z)!=0)
// 			  {
// 			    cout << "warning - imaginary flow" << endl;
// 			    cout << GSL_IMAG(z) << endl;
// 			  }
		      }

		    if(utau<1)
		      {
			utau = 1;
			ux=0;
			uy=0;
		      }

		    if(ux>1 || uy>1)
		      {
			utau = 1;
			ux=0;
			uy=0;
		      }
		 
		    if(eps>0.001)
		      {
			averageuT+=sqrt(ux*ux+uy*uy);
			count ++;
		      }
		    
		    
		    // clean up numerical noise outside the interaction region
 		    if (lat->cells[pos]->getTpA()<1e-10 || lat->cells[pos]->getTpB()<1e-10)
		      {
			eps = 0.;
			utau = 1;
			ux=0;
			uy=0;
		      }

		    lat->cells[pos]->setutau(utau);
		    lat->cells[pos]->setux(ux);
		    lat->cells[pos]->setuy(uy);
		    lat->cells[pos]->setEpsilon(eps);
		    
		    // if(si%6==0 && sj%6==0 && si>0 && sj>0 && param->getWriteOutputs()>0)
		    //   {
		    // 	fout << a*it*dtau << " " << x << " " << y << " " << ux << " " << uy << " " 
		    // 	     << eps << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() << endl;
		    //   }
		    // if(sj == N/2  && param->getWriteOutputs()>0)
		    //   {
		    // 	foutx << a*it*dtau << " " << x << " " << y << " " << ux << " " << uy << " " << "0" << " " 
		    // 	      << eps << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy() << " " << it*dtau*it*dtau*lat->cells[pos]->getTtaueta() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTyy() << " " << it*dtau*it*dtau*lat->cells[pos]->getTetaeta() << endl;
		    // 	if(si == N-1)
		    // 	  {
		    // 	    foutx << endl;
		    // 	  }			    
		    //   }

		  }
	      }
	    
	    gsl_vector_complex_free (eval);
	    gsl_matrix_complex_free (evec);
	    
	    //}
	}
    }


  cout << it*dtau*a << " average u_T=" << averageuT/static_cast<double>(count) << endl;

  //  fout.close();
  //  foutx.close();

  cout << "write outputs: " << param->getWriteOutputs() << endl;
  
  double maxtime = param->getMaxtime(); // maxtime is in fm
  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
  
  double outputTime = 0.2;
  itmax = static_cast<int>(floor(outputTime/(a*dtau)+1e-10));
  cout << "it="<< it << ", itmax=" << itmax << endl;


  // output for hydro
  if(it==itmax && param->getWriteOutputs() > 0)
    {
      double alphas = param->getalphas();
      double g = param->getg();
      double gfactor;
      double Qs;
      int hx = param->getSizeOutput();
      int hy = hx;
      int heta = param->getEtaSizeOutput();
      double hL = param->getLOutput();
      double deta = param->getDetaOutput();
      double c = param->getc();
      double muZero = param->getMuZero();
      double g2mu2A, g2mu2B;
      double PI = param->getPi();

      if(hL>L)
	cout << "WARNING: hydro grid length larger than the computed one." << endl;
      
      int xpos, ypos, xposUp, yposUp, pos1, pos2, pos3;
      double fracx, fracy, x1, x2;
      double xlow, xhigh, ylow, yhigh;
      int pos4;
      double resultE, resultutau, resultux, resultuy;
      double resultT00, resultT0x, resultT0y, resultT0eta;

      double ha;
      ha = hL/static_cast<double>(hx);
      
      stringstream streuH_name;
      streuH_name << "epsilon-u-Hydro" << param->getMPIRank() << ".dat";
      string euH_name;
      euH_name = streuH_name.str();

      ofstream foutEps2(euH_name.c_str(),ios::out); 
      foutEps2 << "# dummy " << 1 << " etamax= " << heta
	       << " xmax= " << hx << " ymax= " << hy << " deta= " << deta 
	       << " dx= " << ha << " dy= " << ha << endl; 
      
      for(int ieta=0; ieta<heta; ieta++) // loop over all positions
	{
	  for(int ix=0; ix<hx; ix++) // loop over all positions
	    {
	      for(int iy=0; iy<hy; iy++)
		{
		  x = -hL/2.+ha*ix;
		  y = -hL/2.+ha*iy;
		  
		  xpos = static_cast<int>(floor((x+L/2.)/a+0.0000000001));
		  ypos = static_cast<int>(floor((y+L/2.)/a+0.0000000001));

		  
		  if(xpos<N-1)
		    xposUp = xpos+1;
		  else
		    xposUp = xpos;
		  
		  if(ypos<N-1)
		    yposUp = ypos+1;
		  else
		    yposUp = ypos;
		  
		  xlow = -L/2.+a*xpos;
		  ylow = -L/2.+a*ypos;
		  
		  xhigh = -L/2.+a*xposUp;
		  yhigh = -L/2.+a*yposUp;
		  
		  fracx = (x-xlow)/ha;
		  
		  pos1 = xpos*N+ypos;
		  pos2 = xposUp*N+ypos;
		  pos3 = xpos*N+yposUp;
		  pos4 = xposUp*N+yposUp;
		  
		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getEpsilon());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getEpsilon());
		  else
		    x2 = 0.;
		  
		  fracy = (y-ylow)/ha;
		  
		  resultE = (1.-fracy)*x1+fracy*x2;
		  
		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getutau())+fracx*(lat->cells[pos2]->getutau());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getutau())+fracx*(lat->cells[pos4]->getutau());
		  else
		    x2 = 0.;
		  
		  resultutau = (1.-fracy)*x1+fracy*x2;

		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getux())+fracx*(lat->cells[pos2]->getux());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getux())+fracx*(lat->cells[pos4]->getux());
		  else
		    x2 = 0.;
		  
		  resultux = (1.-fracy)*x1+fracy*x2;

		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getuy())+fracx*(lat->cells[pos2]->getuy());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getuy())+fracx*(lat->cells[pos4]->getuy());
		  else
		    x2 = 0.;
		  
		  resultuy = (1.-fracy)*x1+fracy*x2;
		  

		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getTtautau())+fracx*(lat->cells[pos2]->getTtautau());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getTtautau())+fracx*(lat->cells[pos4]->getTtautau());
		  else
		    x2 = 0.;
		  
		  resultT00 = (1.-fracy)*x1+fracy*x2;
		  
		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getTtaux())+fracx*(lat->cells[pos2]->getTtaux());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getTtaux())+fracx*(lat->cells[pos4]->getTtaux());
		  else
		    x2 = 0.;
		  
		  resultT0x = (1.-fracy)*x1+fracy*x2;
		  
		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getTtauy())+fracx*(lat->cells[pos2]->getTtauy());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getTtauy())+fracx*(lat->cells[pos4]->getTtauy());
		  else
		    x2 = 0.;
		  
		  resultT0y = (1.-fracy)*x1+fracy*x2;
		  
		  if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
		    x1 = (1-fracx)*(lat->cells[pos1]->getTtaueta())+fracx*(lat->cells[pos2]->getTtaueta());
		  else
		    x1 = 0.;
		  
		  if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
		    x2 = (1-fracx)*(lat->cells[pos3]->getTtaueta())+fracx*(lat->cells[pos4]->getTtaueta());
		  else
		    x2 = 0.;
		  
		  resultT0eta = (1.-fracy)*x1+fracy*x2;
		  
		  
		  if(resultutau<1.)
		    {
		      resultutau=1.;
		      resultux=0.;
		      resultuy=0;
		    }


		  if(param->getRunningCoupling())
		    {
		      
		      // if(pos>0 && pos<(N-1)*N+N-1)
		      // 	{
		      // 	  g2mu2A = lat->cells[pos]->getg2mu2A();
		      // 	}
		      // else 
		      // 	g2mu2A = 0;
		      
		      // if(pos>0 && pos<(N-1)*N+N-1)
		      // 	{
		      // 	  g2mu2B = lat->cells[pos]->getg2mu2B();
		      // 	}
		      // else
		      // 	g2mu2B = 0;
		      
		      // if(param->getRunWithQs()==2)
		      // 	{
		      // 	  if(g2mu2A > g2mu2B)
		      // 	    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      // 	  else
		      // 	    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      // 	}
		      // else if(param->getRunWithQs()==0)
		      // 	{
		      // 	  if(g2mu2A < g2mu2B)
		      // 	    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      // 	  else
		      // 	    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      // 	}
		      // else if(param->getRunWithQs()==1)
		      // 	{
		      // 	  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      // 	}
		      
		      // 3 flavors
		      //	      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		      
		      //gfactor = g*g/(4.*PI*alphas);
		      

		      // run with average Q_s only ! local makes no sense here (stuff has moved in the mean time)
		      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));	  
		      gfactor = g*g/(4.*PI*alphas);
		      //     cout << "Qs=" << Qs << endl;
		      //cout << alphas << endl;
		    }
		  else
		    gfactor = 1.;



		  		  
		  if(abs(0.1973269718*resultE*gfactor) > 0.000000000000001 && abs(x) < L/2. && abs(y) < L/2.)
		    foutEps2 << -(heta-1)/2.*deta+deta*ieta << " " << x << " " << y << " " 
			     << abs(0.1973269718*resultE*gfactor) << " " << resultutau << " " << resultux << " " << resultuy << endl;
		      //		     << " " << 0.1973269718*resultT00*gfactor << " " << 0.1973269718*resultT0x*gfactor
		      //	     << " " << 0.1973269718*resultT0y*gfactor << " " << 0.1973269718*resultT0eta*gfactor << endl; 
		  else
		    foutEps2 << -(heta-1)/2.*deta+deta*ieta << " " << x << " " << y << " " 
			     << 0. << " " << 1. << " " << 0. << " " << 0. << endl;
		  //			     << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl; 
		}
	    }
	  foutEps2 << endl;
	}
      foutEps2.close();
    }
  // done output for hydro
}
