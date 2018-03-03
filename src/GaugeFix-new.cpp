// GaugeFix.cpp is part of the CYM evolution.
// Copyright (C) 2012 Bjoern Schenke.

#include "GaugeFix.h"

//**************************************************************************
// GaugeFix class.

void GaugeFix::FFTChi2(Lattice* lat, Group* group, Parameters *param, int steps)
{
  random->init_genrand64(param->getRandomSeed());
  const int N = param->getSize();
  int nn[2];
  nn[0] = N;
  nn[1] = N;
  int pos, posX, posY, posmX, posmY;
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;
  int maxi=0;
  int maxj=0;
  
  Matrix one(Nc,1.);
  Matrix g(Nc), gdag(Nc), temp(Nc);
  Matrix divA(Nc), divAdag(Nc);
  double energy=0.;
  Matrix Uplaq(Nc), U(Nc), UDx(Nc), UDy(Nc);
  Matrix UmxDag(Nc), UmyDag(Nc);
  double kx, ky, kt2;
  double tr;

  int max_gfiter = steps;
  double maxres = 0.;
  
  Matrix zero(Nc,0.);
  double gresidual = 10000.;
  
  Matrix **chi;
  chi = new Matrix*[N*N];
  
  for(int i=0; i<N*N; i++)
    {
      chi[i] = new Matrix(Nc,0.);
    }
  
  for (int gfiter=0; gfiter<max_gfiter; gfiter++)
    {
      gresidual = 0.;
      maxres = 0.;
      for(int i=0; i<N; i++)
	{
	  for(int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      if(i>0)
		posmX = (i-1)*N+j;
	      else
		posmX = (N-1)*N+j;
	      if(j>0)
		posmY = i*N+(j-1);
	      else
		posmY = i*N+(N-1);
	      
	      divA = (lat->cells[pos]->getUx()-lat->cells[posmX]->getUx()
		      + lat->cells[pos]->getUy()-lat->cells[posmY]->getUy());
 
	      divA = divA - (1./static_cast<double>(Nc))*divA.trace()*one;

	      if(i<N-1)
		posX = (i+1)*N+j;
	      else
		posX = j;
	      if(j<N-1)
		posY = i*N+(j+1);
	      else
		posY = i*N;
	      
	      UDx = lat->cells[posY]->getUx();
	      UDy = lat->cells[pos]->getUy();
	      UDx.conjg();
	      UDy.conjg();
	      
	      Uplaq = lat->cells[pos]->getUx()*(lat->cells[posX]->getUy()*(UDx*UDy));
	   //    lat->cells[pos]->setUplaq(Uplaq);
// 	      if(i==N-1)
// 		lat->cells[pos]->setUplaq(one);

	      energy += 2.*(static_cast<double>(Nc)-(Uplaq.trace()).real());

	//       if(gfiter%10==0 && i==N/2 && j==N/2)
// 		{
// 		  cout << gfiter << ": divA=" << endl << divA << endl << endl;
// 		}
	      
	      *chi[pos] = divA;
	      
	      g = *chi[pos];
	      gdag = *chi[pos];
	      gdag.conjg();
	      
	      gresidual += ((gdag*g).trace()).real()/static_cast<double>(Nc);
	      
	      if(((gdag*g).trace()).real()/static_cast<double>(Nc)>maxres)
		{
		  maxi=i;
		  maxj=j;
		}
	      
	      maxres = max(((gdag*g).trace()).real()/static_cast<double>(Nc),maxres);	  
	      
	      
	      // 	  if (((gdag*g).trace()).real()/static_cast<double>(Nc)>5.)
	      //  	    {
	      //  	      cout << "large residue in cell " << i << " " << j << endl;
	      //  	      cout << "value=" << ((gdag*g).trace()).real()/static_cast<double>(Nc) << endl;
	      //  	    }
	      
	    } // i loop
	} // j loop
      
      gresidual /= (N)*(N);
      energy /= (N)*(N);
      cout << "gauge fixing iteration " << gfiter << endl;
      cout << "residual = " << gresidual << endl;
      cout << "maxres=" << maxres << " at " << maxi << ", " << maxj << endl;
      
      if (gresidual<0.00001*energy)
	break;
      
      fft->fftn(chi,chi,nn,2,1);
      
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      kx = 2.*param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N));
	      ky = 2.*param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N));
	      kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum squared
	      
	      if (kt2!=0.)
	      //	      *chi[pos] = 0.32*(1./(kt2+0.001))*(*chi[pos]); 
	      //*chi[pos] = 0.8*(1./(kt2+0.0000000000001))*(*chi[pos]); 
		*chi[pos] = 0.8*(1./(kt2))*(*chi[pos]); 
	      else 
		*chi[pos] = zero;
	      
	    }
	}
      
      fft->fftn(chi,chi,nn,2,-1);
      
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      
	      g = complex<double>(0.,1.)*(*chi[pos]);
	      
	      // exponentiate
	      g.expm();
	      
	      //	  cout <<  i << " " << j << "  g=" << g.trace() << endl; 
	      
	      // reunitarize
	      g.reu();
	      
	      if(g(2)!=g(2))
		{
		  cout << "problem at " << i << " " << j << " with g=" << g << endl; 
		  g = one;
		}
	      
	      lat->cells[pos]->setg(g);
	      gaugeTransform(lat, group, param, i, j);
	      
	    }
	}
      
  //       for(int i=0; i<N; i++)
  // 	{
  // 	  for(int j=0; j<N; j++)
  // 	    {
  // 	      pos = i*N+j;
  // 	      if(i>0)
  // 		posmX = (i-1)*N+j;
  // 	      else
  // 		posmX = pos;
  // 	      if(j>0)
  // 		posmY = i*N+(j-1);
  // 	      else
  // 		posmY = pos;
  //        //     // reunitarize
  // // 	      U = lat->cells[pos]->getUx();
  // // 	      U.reu();
  // // 	      lat->cells[pos]->setUx(U);
  // // 	      U = lat->cells[pos]->getUy();
  // // 	      U.reu();
  // // 	      lat->cells[pos]->setUy(U);
  // // 	      // check 
  
  // 	      if(i==N/2 && j==N/2)
  // 		{
  // 		  cout << gfiter << ": divA=" << endl << divA << endl << endl;
  // 		  divA = lat->cells[pos]->getUx()-lat->cells[posmX]->getUx()
  // 		    + lat->cells[pos]->getUy()-lat->cells[posmY]->getUy();
  // 		}
  // 	    }
  // 	}
  
} // gfiter loop



for(int i=0; i<N*N; i++)
  {
    delete chi[i];
  }

delete [] chi;
}




void GaugeFix::FFTChi(Lattice* lat, Group* group, Parameters *param, int steps)
{
  const int N = param->getSize();
  int nn[2];
  nn[0] = N;
  nn[1] = N;
  int pos, posX, posY, posmX, posmY;
  int Nc = param->getNc();
  int Nc2m1 = Nc*Nc-1;
  int maxi=0;
  int maxj=0;
  
  Matrix one(Nc,1.);
  Matrix oldg(Nc), g(Nc), gdag(Nc), temp(Nc);
  Matrix divA(Nc), divAdag(Nc);
  Matrix divANew(Nc);
  double energy=0.;
  Matrix Uplaq(Nc), U(Nc), UDx(Nc), UDy(Nc), UDxMx(Nc), UDyMy(Nc) , Ux(Nc), Uy(Nc), UxMx(Nc), UyMy(Nc);
  Matrix UmxDag(Nc), UmyDag(Nc);
  double kx, ky, kt2;
  double tr;

  int max_gfiter = steps;
  double maxres = 0.;
  
  Matrix zero(Nc,0.);
  double gresidual = 10000.;
  
  Matrix **chi;
  chi = new Matrix*[N*N];

  cout << "gauge fixing" << endl;

  for(int i=0; i<N*N; i++)
    {
      chi[i] = new Matrix(Nc,0.);
    }
  
  for (int gfiter=0; gfiter<max_gfiter; gfiter++)
    {
      gresidual = 0.;
      maxres = 0.;
      energy=0.;
      for(int i=0; i<N; i++)
	{
	  for(int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      if(i>0)
		posmX = (i-1)*N+j;
	      else
		posmX = (N-1)*N+j;
	      if(j>0)
		posmY = i*N+(j-1);
	      else
		posmY = i*N+(N-1);
	      
	      if(i<N-1)
		posX = (i+1)*N+j;
	      else
		posX = j;
	      if(j<N-1)
		posY = i*N+(j+1);
	      else
		posY = i*N;
	      
	      Ux = UDx = lat->cells[pos]->getUx();
	      Uy = UDy = lat->cells[pos]->getUy();
	      UxMx = UDxMx = lat->cells[posmX]->getUx();
	      UyMy = UDyMy = lat->cells[posmY]->getUy();
	      UDx.conjg();
	      UDy.conjg();
	      UDxMx.conjg();
	      UDyMy.conjg();
	   
	      //	      divA = (lat->cells[pos]->getUx()-UDx)-(lat->cells[posmX]->getUx()-UDxMx)
	      // 	+(lat->cells[pos]->getUy()-UDy)-(lat->cells[posmY]->getUy()-UDyMy);
	   
	      //divA = Ux-UxMx+UDx-UDxMx+Uy-UyMy+UDy-UDyMy;
	      //divA = Ux-UDx-UxMx+UDxMx+Uy-UDy-UyMy+UDyMy;
	      //divA = (Ux+UDxMx+Uy+UDyMy);
	      divA = (Ux-UxMx+Uy-UyMy);
	      
	      // divANew=zero;
	      // for(int it=0; it<Nc2m1; it++)
	      //  	divANew = divANew+((group->getT(it)*divA).trace())*group->getT(it);

	      g=zero;
	      for (int i=0; i<Nc2m1; i++)
		{
		  g = g + ((divA)*group->getT(i)).trace().imag()*group->getT(i);
		}
	


	      *chi[pos]=gdag=g;// = ( divA - (1./static_cast<double>(Nc))*divA.trace()*one);
	      // divA;
	     
	      // cout << "divA = " << divA << endl << endl;
	      // cout << "divANew = " << divANew << endl << endl;
	      // cout << "tr divA = " << divA.trace() << endl << endl;
	      


// 	      //redefine UDx
// 	      UDx = lat->cells[posY]->getUx();
// 	      UDx.conjg();
	      
// 	      Uplaq = Ux*(lat->cells[posX]->getUy()*(UDx*UDy));
// 	      lat->cells[pos]->setUplaq(Uplaq);
// 	      if(i==N-1)
// 		lat->cells[pos]->setUplaq(one);
	      
// 	      energy += 2.*(static_cast<double>(Nc)-(Uplaq.trace()).real());
		  
	      if(gfiter%100==0 && i==N/2 && j==N/2)
		{
		  // cout << gfiter << ": divA=" << endl << divA << endl << endl;
		  
		  cout << "divA components: " << endl;
		  for (int i=0; i<Nc2m1; i++)
		    {
		      cout << (-1.*(g)*group->getT(i)).trace().real() << " " <<  (-1.*(g)*group->getT(i)).trace().imag()
			   << endl;
		    }
		  cout << endl;
		  for (int i=0; i<Nc2m1; i++)
		    {
		      cout << (-1.*(Ux-UDx-UxMx+UDxMx+Uy-UDy-UyMy+UDyMy)*group->getT(i)).trace().real() << " " <<  (-1.*(Ux-UDx-UxMx+UDxMx+Uy-UDy-UyMy+UDyMy)*group->getT(i)).trace().imag()
			   << endl;
		    }
		  
		}
	      
	      gdag.conjg();
	      
	      gresidual += ((gdag*g).trace()).real()/static_cast<double>(Nc);

	      

// 	      if(((gdag*g).trace()).real()/static_cast<double>(Nc)>maxres)
// 		{
// 		  maxi=i;
// 		  maxj=j;
// 		}
	      
// 	      maxres = max(((gdag*g).trace()).real()/static_cast<double>(Nc),maxres);	  
	      
	      
	      // 	  if (((gdag*g).trace()).real()/static_cast<double>(Nc)>5.)
	      //  	    {
	      //  	      cout << "large residue in cell " << i << " " << j << endl;
	      //  	      cout << "value=" << ((gdag*g).trace()).real()/static_cast<double>(Nc) << endl;
	      //  	    }
	      
	    } // i loop
	} // j loop
      
      gresidual /= (N)*(N);
      //energy /= (N)*(N);
      
      if(gfiter%10==0)
	{
	  //	  cout << "gauge fixing iteration " << gfiter << endl;
	  cout << gfiter << " " << gresidual << endl;
	  //	  cout << "maxres=" << maxres << " at " << maxi << ", " << maxj << endl;
	}
	
      //      if (gresidual<1e-10*energy)
      //      if (gresidual<1e-6*energy)
      if (gresidual<1e-10)
	{
	  break;
	}    

      fft->fftn(chi,chi,nn,2,1);
      
      for (int i=0; i<N; i++)
      	{
      	  for (int j=0; j<N; j++)
      	    {
      	      pos = i*N+j;
      	      kx = sin(param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N)));
      	      ky = sin(param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N)));
      	      kt2 = 4.*(kx*kx+ky*ky); //lattice momentum squared
	      
      	      //     	      if (kt2!=0.)
      	      *chi[pos] = -0.5 * (1./(kt2+1e-9))*(*chi[pos]); 
      	      //  *chi[pos] = - 0.1 * (1./(kt2))*(*chi[pos]); 
      	      //*chi[pos] = (1./(kt2))*(*chi[pos]); 
      	      // else 
      	      //	*chi[pos] = zero;
	      
      	    }
      	}
      
      fft->fftn(chi,chi,nn,2,-1);
      
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      pos = i*N+j;
	      
	      // exponentiate
	      //g = one + complex<double>(1.,0.)*(*chi[pos]);
	      
	      g = complex<double>(0,1.)*(*chi[pos]);
	      g.expm();
	      // cout <<  i << " " << j << "  g=" << g.trace() << endl; 
	      

	      // reunitarize
	      g.reu();

	      if(g(2)!=g(2))
		{
		  cout << "problem at " << i << " " << j << " with g=" << g << endl; 
		  g = one;
		}
	      
	      lat->cells[pos]->setg(g);
	      gaugeTransform(lat, group, param, i, j);
	      
	    }
	}
      
  //       for(int i=0; i<N; i++)
  // 	{
  // 	  for(int j=0; j<N; j++)
  // 	    {
  // 	      pos = i*N+j;
  // 	      if(i>0)
  // 		posmX = (i-1)*N+j;
  // 	      else
  // 		posmX = pos;
  // 	      if(j>0)
  // 		posmY = i*N+(j-1);
  // 	      else
  // 		posmY = pos;
  //        //     // reunitarize
  // // 	      U = lat->cells[pos]->getUx();
  // // 	      U.reu();
  // // 	      lat->cells[pos]->setUx(U);
  // // 	      U = lat->cells[pos]->getUy();
  // // 	      U.reu();
  // // 	      lat->cells[pos]->setUy(U);
  // // 	      // check 
  
  // 	      if(i==N/2 && j==N/2)
  // 		{
  // 		  cout << gfiter << ": divA=" << endl << divA << endl << endl;
  // 		  divA = lat->cells[pos]->getUx()-lat->cells[posmX]->getUx()
  // 		    + lat->cells[pos]->getUy()-lat->cells[posmY]->getUy();
  // 		}
  // 	    }
  // 	}
  
} // gfiter loop



for(int i=0; i<N*N; i++)
  {
    delete chi[i];
  }

delete [] chi;
}

	


void GaugeFix::gaugeTransform(Lattice* lat, Group* group, Parameters *param, int i, int j)
{
  int N = param->getSize();
  int pos, posmX, posmY;
  int Nc = param->getNc();
  Matrix g(Nc), gdag(Nc), temp(Nc);

  pos = i*N+j;
  if(i>0)
    posmX = (i-1)*N+j;
  else
    posmX = (N-1)*N+j;

  if(j>0)
    posmY = i*N+(j-1);
  else 
    posmY = i*N+N-1;

  g = gdag = lat->cells[pos]->getg();
  gdag.conjg();


//   if (i==N/2 && j==N/2)
//     cout << "g=" << g << endl;
  
      // gauge transform Ux and Uy
      lat->cells[pos]->setUx( g * lat->cells[pos]->getUx() );
      lat->cells[pos]->setUy( g * lat->cells[pos]->getUy() );
      
      lat->cells[posmX]->setUx( lat->cells[posmX]->getUx() * gdag );
      lat->cells[posmY]->setUy( lat->cells[posmY]->getUy() * gdag );

      // gauge transform Ex and Ey
      lat->cells[pos]->setE1( g * lat->cells[pos]->getE1() * gdag );
      lat->cells[pos]->setE2( g * lat->cells[pos]->getE2() * gdag );
      
      // gauge transform phi and pi
      lat->cells[pos]->setphi( g * lat->cells[pos]->getphi() * gdag  );
      lat->cells[pos]->setpi( g * lat->cells[pos]->getpi() * gdag );
}
