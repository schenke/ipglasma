// Evolution.cpp is part of the IP-Glasma evolution solver.
// Copyright (C) 2012 Bjoern Schenke.
#include "Evolution.h"

//**************************************************************************
// Evolution class.

void Evolution::evolveU(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau)
{
  // tau is the current time. The time argument of E^i is tau+dtau/2
  // we evolve to tau+dtau
  const int Nc = param->getNc();
  const int N = param->getSize();
  const double g = param->getg();

  const int n = 2;
  const Matrix one(Nc,1.);

#pragma omp parallel
  {
    Matrix E1(Nc);
    Matrix E2(Nc);
    
    Matrix temp1(Nc);
    Matrix temp2(Nc);

#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        // retrieve current E1 and E2 (that's the one defined at half a time step in the future (from tau))
        E1 = complex<double>(0.,g*g*dtau/(tau+dtau/2.))*lat->cells[pos]->getE1();
        //E1.expm(); // E1 now contains the exponential of i g^2 dtau/(tau+dtau/2)*E1
        
        temp2 = one + 1./(double)n * E1;
        for (int in=0; in<n-1; in++) 
          {
            temp1 = E1*temp2;
            temp2 = one + 1./(double)(n-1-in) * temp1;
          }
        
        E1 = temp2;
        
        E2 = complex<double>(0.,g*g*dtau/(tau+dtau/2.))*lat->cells[pos]->getE2();
        // E2.expm(); // E2 now contains the exponential of i g^2 dtau/(tau+dtau/2)*E2
        
        temp2 = one + 1./(double)n * E2;
        for (int in=0; in<n-1; in++) 
          {
            temp1 = E2*temp2;
            temp2 = one + 1./(double)(n-1-in) * temp1;
          }
        
        E2 = temp2;
        
        bufferlat->cells[pos]->setbuffer1(E1*lat->cells[pos]->getUx()); 
        bufferlat->cells[pos]->setbuffer2(E2*lat->cells[pos]->getUy()); 
      }
  
 #pragma omp for  
    for (int pos=0; pos<N*N; pos++)
      {
        lat->cells[pos]->setUx(bufferlat->cells[pos]->getbuffer1()); 
        lat->cells[pos]->setUy(bufferlat->cells[pos]->getbuffer2()); 
      }
  } 
}
  

void Evolution::evolvePhi(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau)
{
  // tau is the current time. The time argument of pi is tau+dtau/2
  // we evolve to tau+dtau
  const int Nc = param->getNc();
  const int N = param->getSize();
  const double g = param->getg();


#pragma omp parallel
  {
    Matrix phi(Nc);
    Matrix pi(Nc);
    
#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        // retrieve current phi (at time tau)
        phi = lat->cells[pos]->getphi();
        // retrieve current pi (at time tau+dtau/2)
        pi = lat->cells[pos]->getpi();
        
        phi = phi + (tau+dtau/2.)*dtau*pi;
        
        //set the new phi (at time tau+dtau)
        bufferlat->cells[pos]->setbuffer1(phi); 
      }
 
#pragma omp for   
    for (int pos=0; pos<N*N; pos++)
      {
        lat->cells[pos]->setphi(bufferlat->cells[pos]->getbuffer1()); 
      }
  }
}

void Evolution::evolvePi(Lattice* lat, BufferLattice * bufferlat, Group* group, Parameters *param, double dtau, double tau)
{
  const int Nc = param->getNc();
  const int N = param->getSize();
  const double g = param->getg();

#pragma omp parallel 
  {
    Matrix Ux(Nc);
    Matrix Uy(Nc);
    Matrix UxXm1(Nc);
    Matrix UyYm1(Nc);
    
    Matrix phi(Nc);
    Matrix phiX(Nc); // this is \tilde{phi}_x
    Matrix phiY(Nc); // this is \tilde{phi}_y
    Matrix phimX(Nc); // this is \tilde{-phi}_x
    Matrix phimY(Nc); // this is \tilde{-phi}_y
    Matrix bracket(Nc); // this will hold [phiX+phimX-2*phi+phiY+phimY-2*phi]
    Matrix pi(Nc);
    Matrix one(Nc,1.);  
    // Matrix zero(Nc,0.);  
    
#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        // retrieve current Ux and Uy and compute conjugates
        Ux = lat->cells[pos]->getUx();
        Uy = lat->cells[pos]->getUy();
        // retrieve current pi (at time tau-dtau/2)
        pi = lat->cells[pos]->getpi();
        // retrieve current phi (at time tau) at this x_T
        phi = lat->cells[pos]->getphi();
        
        // retrieve current phi (at time tau) at x_T+1
        // parallel transport:
        phiX = Ux*Ux.prodABconj(lat->cells[lat->pospX[pos]]->getphi(),Ux);
        phiY = Uy*Uy.prodABconj(lat->cells[lat->pospY[pos]]->getphi(),Uy);
        
        // phi_{-x} should be defined as UxD*phimX*Ux with the Ux and UxD reversed from the phi_{+x} case
        // retrieve current phi (at time tau) at x_T-1
        // parallel transport:
        UxXm1 = lat->cells[lat->posmX[pos]]->getUx();
        UyYm1 = lat->cells[lat->posmY[pos]]->getUy();
        
        phimX = Ux.prodAconjB(UxXm1,lat->cells[lat->posmX[pos]]->getphi())*UxXm1;
        phimY = Ux.prodAconjB(UyYm1,lat->cells[lat->posmY[pos]]->getphi())*UyYm1;
        
        bracket = phiX + phimX + phiY + phimY - 4.*phi; // sum over both directions is included here 
        
        //	  if (i>0 && i<N-1 && j>0 && j<N-1)
        pi += dtau/(tau)*bracket; // divide by \tau because this is computing pi(tau+dtau/2) from pi(tau-dtau/2) and phi(tau)
        
        //pi = pi - pi.trace()/static_cast<double>(Nc)*one;
        
        // set the new pi (at time tau+dtau/2)

        bufferlat->cells[pos]->setbuffer1(pi); 

      }
#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        lat->cells[pos]->setpi(bufferlat->cells[pos]->getbuffer1()); 
      }
  }
}

void Evolution::evolveE(Lattice* lat, BufferLattice *bufferlat, Group* group, Parameters *param, double dtau, double tau)
{
  const int Nc = param->getNc();
  const int N = param->getSize();
  const double g = param->getg();

#pragma omp parallel
  {
    Matrix Ux(Nc);
    Matrix Uy(Nc);
    Matrix temp1(Nc); // can contain p or m
    Matrix temp2(Nc);
    Matrix temp3(Nc);
    Matrix En(Nc);
    Matrix phi(Nc);
    Matrix phiN(Nc); // this is \tilde{phi}_x OR \tilde{phi}_y
    
    //plaquettes:
    Matrix U12(Nc);
    Matrix U1m2(Nc);
    Matrix U12Dag(Nc); // equals (U21)
    Matrix U2m1(Nc);
    complex<double> trace;
    Matrix one(Nc,1.);
    
#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
        En = lat->cells[pos]->getE1();
        // retrieve current phi (at time tau) at this x_T
        phi = lat->cells[pos]->getphi();
        // retrieve current phi (at time tau) at x_T+1
        phiN = lat->cells[lat->pospX[pos]]->getphi();
        // parallel transport:
        // retrieve current Ux and Uy
        Ux = lat->cells[pos]->getUx();
        phiN = Ux*Ux.prodABconj(phiN,Ux);
        
        // compute plaquettes:
        Uy = lat->cells[pos]->getUy();	  
        temp1 = lat->cells[lat->pospY[pos]]->getUx(); //UxYp1Dag
        temp1.conjg();
        U12 = (Ux*lat->cells[lat->pospX[pos]]->getUy())*(Ux.prodABconj(temp1,Uy));
        
        temp1 = lat->cells[lat->posmY[pos]]->getUx(); //UxYm1Dag
        temp2 = lat->cells[lat->pospXmY[pos]]->getUy(); //UyXp1Ym1Dag
        U1m2 = (Ux.prodABconj(Ux,temp2))*(Ux.prodAconjB(temp1,lat->cells[lat->posmY[pos]]->getUy()));
        
        temp1 = lat->cells[lat->posmX[pos]]->getUy(); //UyXm1Dag
        temp2 = lat->cells[lat->posmXpY[pos]]->getUx(); //UxXm1Yp1Dag
        U2m1 = (Ux.prodABconj(Uy,temp2))*(Ux.prodAconjB(temp1,lat->cells[lat->posmX[pos]]->getUx()));
        
        U12Dag = U12;
        U12Dag.conjg();
        
        // do E1 update:
        
        temp3 = U1m2;
        temp3.conjg();
        
        temp1 = U12;
        temp1 += U1m2;
        temp1 -= U12Dag;
        temp1 -= temp3;
        
        trace = temp1.trace();
        temp1 -= trace/static_cast<double>(Nc)*one;
        
        temp2 = phiN*phi - phi*phiN;
        
        //	  if (i>0 && i<N-1 && j>0 && j<N-1)
        En += complex<double>(0.,1.)*tau*dtau/(2.*g*g)*temp1 + complex<double>(0.,1.)*dtau/tau*temp2;
        
        trace = En.trace();
        En -= (trace/static_cast<double>(Nc))*one;
        bufferlat->cells[pos]->setbuffer1(En); 
       
        
        // do E2 update:
        
        temp3 = U2m1;
        temp3.conjg();
        
        temp1 = U12Dag;
        temp1 += U2m1;
        temp1 -= U12;
        temp1 -= temp3;
        trace = temp1.trace();
        temp1 -= (trace/static_cast<double>(Nc))*one;
        
        phiN = lat->cells[lat->pospY[pos]]->getphi();
        phiN = Uy*Uy.prodABconj(phiN,Uy);
        
        temp2 = phiN*phi - phi*phiN;
        
        En = lat->cells[pos]->getE2();
        //if (i>0 && i<N-1 && j>0 && j<N-1)
        En += complex<double>(0.,1.)*tau*dtau/(2.*g*g)*temp1 + complex<double>(0.,1.)*dtau/tau*temp2;
        
        trace = En.trace();
        En -= (trace/static_cast<double>(Nc))*one;
        bufferlat->cells[pos]->setbuffer2(En); 
      }

#pragma omp for
    for (int pos=0; pos<N*N; pos++)
      {
        lat->cells[pos]->setE1(bufferlat->cells[pos]->getbuffer1());
        lat->cells[pos]->setE2(bufferlat->cells[pos]->getbuffer2());
      }
  }
}

void Evolution::checkGaussLaw(Lattice* lat, Group* group, Parameters *param, double dtau, double tau)
{
  const int Nc = param->getNc();
  const int N = param->getSize();
  const double g = param->getg();

  Matrix Ux(Nc);
  Matrix UxXm1(Nc);
  Matrix UxYp1(Nc);
  Matrix UxYm1(Nc);
  Matrix UxXm1Yp1(Nc);
  
  Matrix Uy(Nc);
  Matrix UyXm1(Nc);
  Matrix UyXp1(Nc);
  Matrix UyYm1(Nc);
  Matrix UyXp1Ym1(Nc);
  
  Matrix UxDag(Nc);
  Matrix UxXm1Dag(Nc);
  Matrix UxYp1Dag(Nc);
  Matrix UxYm1Dag(Nc);
  Matrix UxXm1Yp1Dag(Nc);
  
  Matrix UyDag(Nc);
  Matrix UyXm1Dag(Nc);
  Matrix UyXp1Dag(Nc);
  Matrix UyYm1Dag(Nc);
  Matrix UyXp1Ym1Dag(Nc);
  
  Matrix E1(Nc);
  Matrix E2(Nc);
  Matrix E1mX(Nc);
  Matrix E2mY(Nc);
  Matrix phi(Nc);
  Matrix pi(Nc);
  Matrix phiX(Nc); // this is \tilde{phi}_x
  Matrix phiY(Nc); // this is \tilde{phi}_y
  
  Matrix Gauss(Nc);
  complex<double> trace;
  Matrix one(Nc,1.);
  double largest=0;
  
  for (int pos=0; pos<N*N; pos++)
    {
      // retrieve current Ux and Uy
      Ux = lat->cells[pos]->getUx();
      Uy = lat->cells[pos]->getUy();
      UxDag = Ux;
      UxDag.conjg();
      UyDag = Uy;
      UyDag.conjg();
      
      UxXm1 = lat->cells[lat->posmX[pos]]->getUx();
      UxYm1 = lat->cells[lat->posmY[pos]]->getUx();
      UxXm1Dag = UxXm1;
      UxXm1Dag.conjg();
      UxYm1Dag = UxYm1;
      UxYm1Dag.conjg();
      
      UyXm1 = lat->cells[lat->posmX[pos]]->getUy();
      UyYm1 = lat->cells[lat->posmY[pos]]->getUy();
      UyXm1Dag = UyXm1;
      UyXm1Dag.conjg();
      UyYm1Dag = UyYm1;
      UyYm1Dag.conjg();
      
      // retrieve current E1 and E2 (that's the one defined at tau-dtau/2)
      E1 = lat->cells[pos]->getE1();
      E2 = lat->cells[pos]->getE2();
      E1mX = lat->cells[lat->posmX[pos]]->getE1();
      E2mY = lat->cells[lat->posmY[pos]]->getE2();
      // retrieve current phi (at time tau) at this x_T
      phi = lat->cells[pos]->getphi();
      // retrieve current pi
      pi = lat->cells[pos]->getpi();
      
      Gauss = UxXm1Dag*E1mX*UxXm1 - E1 + UyYm1Dag*E2mY*UyYm1 - E2 -complex<double>(0.,1.)*(phi*pi-pi*phi);
      
      if (Gauss.square()>largest)
        largest = Gauss.square();
      
    }
  cout << "Gauss violation=" << largest << endl;
}

void Evolution::run(Lattice* lat, BufferLattice* bufferlat, Group* group, Parameters *param)
{
  int nn[2];
  int Nc = param->getNc();
  int pos, pos1, pos2, pos3, posx, posy, posxm, posym, posxmym;
  int posX, posY, posXY;
  int counts, countMe;
  int N = param->getSize();
  nn[0] = N;
  nn[1] = N;
  double g = param->getg();
  double epsilon;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double avgEps;
  double avgBl2;
  double avgEl2;
  double avgEt2;
  double avgBt2;
  double x,y;

  Matrix one(Nc,1.);

  Matrix zero(Nc,0.);
  Matrix AM(int(Nc),0.);
  Matrix AP(int(Nc),0.);
  Matrix Uplaq(int(Nc),0.);
  Matrix Uplaq1(int(Nc),0.);
  Matrix Uplaq2(int(Nc),0.);
  Matrix Uplaq3(int(Nc),0.);
  Matrix Uplaq4(int(Nc),0.);
  Matrix Ux(int(Nc),0.);
  Matrix Uy(int(Nc),0.);
  Matrix UxDag(int(Nc),0.);
  Matrix UyDag(int(Nc),0.);
  Matrix Ux2(int(Nc),0.);
  Matrix Uy2(int(Nc),0.);
  Matrix UD(int(Nc),0.);
  Matrix UDx(int(Nc),0.);
  Matrix UDy(int(Nc),0.);
  Matrix UDx1(int(Nc),0.);
  Matrix UDy1(int(Nc),0.);
  Matrix E1(int(Nc),0.);
  Matrix E2(int(Nc),0.);
  Matrix phi(int(Nc),0.);
  Matrix phiX(int(Nc),0.);
  Matrix phiY(int(Nc),0.);
  Matrix pi(Nc);

  double alphas;
  double gfactor;
  double Qs, g2mu2A, g2mu2B;
  double muZero = param->getMuZero();
  double c = param->getc();
  Matrix phiTildeX(Nc);
  Matrix phiTildeY(Nc);

  
  // stringstream stru_name;
  // stru_name << "u" << param->getMPIRank() << ".dat";
  // string u_name;
  // u_name = stru_name.str();
  // ofstream foutu(u_name.c_str(),ios::out); 
  // foutu.close();

  // stringstream stru3D_name;
  // stru3D_name << "u3D" << param->getMPIRank() << ".dat";
  // string u3D_name;
  // u3D_name = stru3D_name.str();
  // ofstream foutu3D(u3D_name.c_str(),ios::out); 
  // foutu3D.close();

  // stringstream streEP_name;
  // streEP_name << "epsilonEvolutionPlot" << param->getMPIRank() << ".dat";
  // string eEP_name;
  // eEP_name = streEP_name.str();
  // ofstream foutee(eEP_name.c_str(),ios::out); 
  // foutee.close();

  // stringstream strux_name;
  // strux_name << "u-x" << param->getMPIRank() << ".dat";
  // string ux_name;
  // ux_name = strux_name.str();
  // ofstream foutux(ux_name.c_str(),ios::out); 
  // foutux.close();

  // stringstream strux3D_name;
  // strux3D_name << "u-x3D" << param->getMPIRank() << ".dat";
  // string ux3D_name;
  // ux3D_name = strux3D_name.str();
  // ofstream foutux3D(ux3D_name.c_str(),ios::out); 
  // foutux3D.close();

  // stringstream strNpartdEdy_name;
  // strNpartdEdy_name << "NpartdEdy" << param->getMPIRank() << ".dat";
  // string NpartdEdy_name;
  // NpartdEdy_name = strNpartdEdy_name.str();

  // stringstream strdEdy_name;
  // strdEdy_name << "dEdy" << param->getMPIRank() << ".dat";
  // string dEdy_name;
  // dEdy_name = strdEdy_name.str();

  // stringstream strEB_name;
  // strEB_name << "E2B2" << param->getMPIRank() << ".dat";
  // string EB_name;
  // EB_name = strEB_name.str();
  // ofstream foutEB(EB_name.c_str(),ios::out); 
   

  // do the first half step of the momenta (E1,E2,pi)
  // for now I use the \tau=0 value at \tau=d\tau/2.

  double dtau = param->getdtau(); // dtau is in lattice units
  
  double maxtime;
  if ( param->getInverseQsForMaxTime() == 1 )
    {
      maxtime = 1./param->getAverageQs()*hbarc;
      cout << "maximal evolution time = " << maxtime << " fm" << endl; 
    }
  else
    {
      maxtime = param->getMaxtime(); // maxtime is in fm
    }

  // stringstream strame_name;
  // strame_name << "AverageMaximalEpsilon" << param->getMPIRank() << ".dat";
  // string ame_name;
  // ame_name = strame_name.str();

  // ofstream foutEpsA(ame_name.c_str(),ios::app); 


  //  cout << "check" << endl;

  // E and Pi at tau=dtau/2 are equal to the initial ones (at tau=0)
  // now evolve phi and U to time tau=dtau.
  evolvePhi(lat, bufferlat, group, param, dtau, 0.);
  evolveU(lat, bufferlat, group, param, dtau, 0.);

  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
  int it1 = static_cast<int>(floor(0.2/(a*dtau)+1e-10));
  int it2 = static_cast<int>(floor(0.4/(a*dtau)+1e-10));
    
  
  //  Tmunu(lat,group,param,0);

  cout << "Starting evolution" << endl;
  cout << "itmax=" << itmax << endl;


  // do evolution
  for (int it=1; it<=itmax; it++)
    {
      if(it%10==0)
	cout << endl << "Evolving to time " << it*a*dtau << " fm/c" << endl << endl;  
      // evolve from time tau-dtau/2 to tau+dtau/2 
      if (it<itmax)
	{
	  evolvePi(lat, bufferlat, group, param, dtau, (it)*dtau); // the last argument is the current time tau. 
	  evolveE(lat, bufferlat, group, param, dtau, (it)*dtau);
	  
	  // evolve from time tau to tau+dtau
	  evolvePhi(lat, bufferlat, group, param, dtau, (it)*dtau);
	  evolveU(lat, bufferlat, group, param, dtau, (it)*dtau);
	}
      else if(it==itmax)
	{
	  evolvePi(lat, bufferlat, group, param, dtau/2., (it)*dtau); // the last argument is the current time tau. 
	  evolveE(lat, bufferlat, group, param, dtau/2., (it)*dtau);
  	}

      //      Tmunu(lat,group,param,it);


      if(it == 1 || it==floor(it1) || it==floor(it2) || it==floor(itmax))
	{	  
	  Tmunu(lat,group,param,it);
	  //eccentricity(lat,group,param,it,1.,1);
	  if (it > 1)
	    u(lat,group,param,it); // computes flow velocity and correct energy density 
	}

      // avgEps=0.;
      // avgBt2=0.;
      // avgBl2=0.;
      // avgEt2=0.;
      // avgEl2=0.;
      // double dN = static_cast<double>(N);
      
      // countMe = 0;
      // //cout << "output energydensity" << endl;
      // // output epsilon 

      // // double maxtime = param->getMaxtime(); // maxtime is in fm
      // //int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
      
      

      // // output for plots
      // double g2mu = param->getg2mu();
      // for(int i=0; i<N; i++) // loop over all positions
      // 	{
      // 	  for(int j=0; j<N; j++)
      // 	    {
	      
      // 	      pos = i*N+j;
      // 	      if(i<N-1)
      // 		posX = (i+1)*N+j;
      // 	      else
      // 		posX = j;
      // 	      if(j<N-1)
      // 		posY = i*N+(j+1);
      // 	      else
      // 		posY = i*N;
	      
      // 	      if(i<N-1 && j<N-1)
      // 		posXY =  (i+1)*N+j+1;
      // 	      else if(i<N-1)
      // 		posXY =  (i+1)*N;
      // 	      else if(j<N-1)
      // 		posXY = j+1;
      // 	      else
      // 		posXY = 0;
	      
      // 	      UDx = lat->cells[posY]->getUx();
      // 	      UDy = lat->cells[pos]->getUy();
      // 	      UDx.conjg();
      // 	      UDy.conjg();
	      
      // 	      Uplaq = lat->cells[pos]->getUx()*(lat->cells[posX]->getUy()*(UDx*UDy));
      // 	      lat->cells[pos]->setUplaq(Uplaq);
      // 	      if(i==N-1)
      // 		lat->cells[pos]->setUplaq(one);

      // 	      phi = lat->cells[pos]->getphi();
	      
      // 	      Ux = lat->cells[pos]->getUx();
      // 	      Uy = lat->cells[pos]->getUy();
      // 	      UDx = Ux;
      // 	      UDx.conjg();
      // 	      UDy = Uy;
      // 	      UDy.conjg();
	      
      // 	      phiX = lat->cells[posX]->getphi();
      // 	      phiY = lat->cells[posY]->getphi();

      // 	      phiTildeX = Ux*phiX*UDx;
      // 	      phiTildeY = Uy*phiY*UDy;      
	      
      // 	      E1 = lat->cells[pos]->getE1();
      // 	      E2 = lat->cells[pos]->getE2();
      // 	      pi = lat->cells[pos]->getpi();  
	      
      // 	      avgEt2 += 1./(it*dtau)/(it*dtau)*real((E1*E1).trace()+(E2*E2).trace());
      // 	      avgEl2 += ((pi*pi).trace()).real();
      // 	      avgBt2 += 1./(it*dtau)/(it*dtau)*(real(((phi-phiTildeX)*(phi-phiTildeX)).trace())
      // 						+real(((phi-phiTildeY)*(phi-phiTildeY)).trace()));
      // 	      avgBl2 += 2*(static_cast<double>(Nc)-(Uplaq.trace()).real());


      // 	    }
      // 	}

      // avgEt2/=(dN*dN);
      // avgEl2/=(dN*dN);
      // avgBt2/=(dN*dN);
      // avgBl2/=(dN*dN);

      

      // // write B^2 and E^2 output here...
      
      
      // //      foutEB << it*dtau << " " << dtau*it*avgEt2 << " " << dtau*it*avgEl2 << " " << dtau*it*avgBt2 << " " << dtau*it*avgBl2 << endl;
    

      if(it==1 && param->getWriteOutputs() == 2)
	{	  
	  stringstream streI_name;
	  streI_name << "epsilonInitialPlot" << param->getMPIRank() << ".dat";
	  string eI_name;
	  eI_name = streI_name.str();

	  ofstream foutEps(eI_name.c_str(),ios::out); 
	  for(int ix=0; ix<N; ix++) // loop over all positions
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = ix*N+iy;
		  x = -L/2.+a*ix;
		  y = -L/2.+a*iy;

		  if(param->getRunningCoupling())
		    {
		      if(pos>0 && pos<(N-1)*N+N-1)
			{
			  g2mu2A = lat->cells[pos]->getg2mu2A();
			}
		      else 
			g2mu2A = 0;
		      
		      if(pos>0 && pos<(N-1)*N+N-1)
			{
			  g2mu2B = lat->cells[pos]->getg2mu2B();
			}
		      else
			g2mu2B = 0;
		      
		      if(param->getRunWithQs()==2)
			{
			  if(g2mu2A > g2mu2B)
			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  else
			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}
		      else if(param->getRunWithQs()==0)
			{
			  if(g2mu2A < g2mu2B)
			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  else
			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}
		      else if(param->getRunWithQs()==1)
			{
			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}
		      
		      // 3 flavors
		      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		      
		      gfactor = g*g/(4.*PI*alphas);
		      //     cout << "Qs=" << Qs << endl;
		      //cout << alphas << endl;
		      // run with the local (in transverse plane) coupling
		    }
		  else
		    gfactor = 1.;

		  foutEps << x << " " << y << " " << 0.1973269718*gfactor*abs(lat->cells[pos]->getEpsilon()) << endl; 
		  // abs just to get rid of negative 10^(-17) numbers at edge
		}
	      foutEps << endl;
	    }
	  foutEps.close();
 	}


      if(it==itmax/2 && param->getWriteOutputs() == 2)
	{	  
	  stringstream streInt_name;
	  streInt_name << "epsilonIntermediatePlot" << param->getMPIRank() << ".dat";
	  string eInt_name;
	  eInt_name = streInt_name.str();

	  ofstream foutEps2(eInt_name.c_str(),ios::out); 
	  for(int ix=0; ix<N; ix++) // loop over all positions
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = ix*N+iy;
		  x = -L/2.+a*ix;
		  y = -L/2.+a*iy;

		  if(param->getRunningCoupling())
		    {
		      if(pos>0 && pos<(N-1)*N+N-1)
			{
			  g2mu2A = lat->cells[pos]->getg2mu2A();
			}
		      else 
			g2mu2A = 0;
		      
		      if(pos>0 && pos<(N-1)*N+N-1)
			{
			  g2mu2B = lat->cells[pos]->getg2mu2B();
			}
		      else
			g2mu2B = 0;
		      
		      if(param->getRunWithQs()==2)
			{
			  if(g2mu2A > g2mu2B)
			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  else
			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}
		      else if(param->getRunWithQs()==0)
			{
			  if(g2mu2A < g2mu2B)
			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			  else
			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}
		      else if(param->getRunWithQs()==1)
			{
			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
			}

		      if (param->getRunWithLocalQs() == 1)
			{
			  // 3 flavors
			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
			  gfactor = g*g/(4.*PI*alphas);
			  //     cout << "Qs=" << Qs << endl;
			  //cout << alphas << endl;
			  // run with the local (in transverse plane) coupling
			}
		      else
 			{
			  if ( param->getRunWithQs() == 0 )
			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
			  else if ( param->getRunWithQs() == 1 )
			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
			  else if ( param->getRunWithQs() == 2 )
			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));

			  gfactor = g*g/(4.*PI*alphas);
			}
		    }
		  else
		    gfactor = 1.;

		  foutEps2 << x << " " << y << " " << 0.1973269718*gfactor*abs(lat->cells[pos]->getEpsilon()) << endl; 
		  // abs just to get rid of negative 10^(-17) numbers at edge
		}
	      foutEps2 << endl;
	    }
	  foutEps2.close();
 	}

    //   if(it==itmax)
    // 	{	  
    // 	  stringstream strePl_name;
    // 	  strePl_name << "epsilonPlot" << param->getMPIRank() << ".dat";
    // 	  string ePl_name;
    // 	  ePl_name = strePl_name.str();

    // 	  ofstream foutEps3(ePl_name.c_str(),ios::app); 
    // 	  for(int ix=0; ix<N; ix++) // loop over all positions
    // 	    {
    // 	      for(int iy=0; iy<N; iy++)
    // 		{
    // 		  pos = ix*N+iy;
    // 		  x = -L/2.+a*ix;
    // 		  y = -L/2.+a*iy;

    // 		  if(param->getRunningCoupling())
    // 		    {
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2A = lat->cells[pos]->getg2mu2A();
    // 			}
    // 		      else 
    // 			g2mu2A = 0;
		      
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2B = lat->cells[pos]->getg2mu2B();
    // 			}
    // 		      else
    // 			g2mu2B = 0;
		      
    // 		      if(param->getRunWithQs()==2)
    // 			{
    // 			  if(g2mu2A > g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==0)
    // 			{
    // 			  if(g2mu2A < g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==1)
    // 			{
    // 			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
		      

    // 		      if ( param->getRunWithLocalQs() == 1) 
    // 			{
    // 			  // 3 flavors
    // 			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
    // 			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			  //     cout << "Qs=" << Qs << endl;
    // 			  //cout << alphas << endl;
    // 			  // run with the local (in transverse plane) coupling
    // 			}
    // 		      else
    // 			{
    // 			  if ( param->getRunWithQs() == 0 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 1 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 2 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
			  
    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			}
    // 		    }
    //              else
    //            gfactor = 1.;
  

    //            if(param->getWriteOutputs() == 2)
    //              foutEps3 << x << " " << y << " " << 0.1973269718*gfactor*abs(lat->cells[pos]->getEpsilon()) << endl;   
    //              // abs just to get rid of negative 10^(-17) numbers at edge
    // 		}
    // 	      foutEps3 << endl;
    // 	    }
    // 	  foutEps3.close();
    // 	}




    // if(it==itmax)
    // 	{	  
    // 	  stringstream strepsy_name;
    // 	  strepsy_name << "eps-y" << param->getMPIRank() << ".dat";
    // 	  string epsy_name;
    // 	  epsy_name = strepsy_name.str();

    // 	  ofstream foutEpsY(epsy_name.c_str(),ios::out); 
    // 	  for(int ix=N/2; ix<=N/2; ix++) // loop over all positions
    // 	    {
    // 	      for(int iy=0; iy<N; iy++)
    // 		{
    // 		  pos = ix*N+iy;
    // 		  x = -L/2.+a*ix;
    // 		  y = -L/2.+a*iy;
		  
    // 		  if(param->getRunningCoupling())
    // 		    {
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2A = lat->cells[pos]->getg2mu2A();
    // 			}
    // 		      else 
    // 			g2mu2A = 0;
		      
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2B = lat->cells[pos]->getg2mu2B();
    // 			}
    // 		      else
    // 			g2mu2B = 0;
		      
    // 		      if(param->getRunWithQs()==2)
    // 			{
    // 			  if(g2mu2A > g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==0)
    // 			{
    // 			  if(g2mu2A < g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==1)
    // 			{
    // 			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
		      

    // 		      if ( param->getRunWithLocalQs() == 1) 
    // 			{
    // 			  // 3 flavors
    // 			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
    // 			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			  //     cout << "Qs=" << Qs << endl;
    // 			  //cout << alphas << endl;
    // 			  // run with the local (in transverse plane) coupling
    // 			}
    // 		      else
    // 			{
    // 			  if ( param->getRunWithQs() == 0 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 1 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 2 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
			  
    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			}
    // 		    }
    //              else
    // 		   gfactor = 1.;
  
    // 		  foutEpsY << x << " " << y << " " << 0.1973269718*gfactor*abs(lat->cells[pos]->getEpsilon()) << endl;   
    //              // abs just to get rid of negative 10^(-17) numbers at edge
    // 		}
    // 	    }
    // 	  foutEpsY.close();
    // 	}



    // if((it==1 || (it%10==0 && it>0))  && param->getWriteOutputs() == 2)
    // 	{	  
    // 	  ofstream foutEps4(eEP_name.c_str(),ios::app); 
    // 	  for(int ix=0; ix<N; ix++) // loop over all positions
    // 	    {
    // 	      for(int iy=0; iy<N; iy++)
    // 		{
    // 		  pos = ix*N+iy;
    // 		  x = -L/2.+a*ix;
    // 		  y = -L/2.+a*iy;

    // 		  if(param->getRunningCoupling())
    // 		    {
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2A = lat->cells[pos]->getg2mu2A();
    // 			}
    // 		      else 
    // 			g2mu2A = 0;
		      
    // 		      if(pos>0 && pos<(N-1)*N+N-1)
    // 			{
    // 			  g2mu2B = lat->cells[pos]->getg2mu2B();
    // 			}
    // 		      else
    // 			g2mu2B = 0;
		      
    // 		      if(param->getRunWithQs()==2)
    // 			{
    // 			  if(g2mu2A > g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==0)
    // 			{
    // 			  if(g2mu2A < g2mu2B)
    // 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			  else
    // 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
    // 		      else if(param->getRunWithQs()==1)
    // 			{
    // 			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
    // 			}
		      
    // 		      if ( param->getRunWithLocalQs() == 1 )
    // 			{
    // 			  // 3 flavors
    // 			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
    // 			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			  //     cout << "Qs=" << Qs << endl;
    // 			  //cout << alphas << endl;
    // 			  // run with the local (in transverse plane) coupling
    // 			}
    // 		      else
    // 			{
    // 			  if ( param->getRunWithQs() == 0 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 1 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
    // 			  else if ( param->getRunWithQs() == 2 )
    // 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));

    // 			  gfactor = g*g/(4.*PI*alphas);
    // 			}
    // 		    }
    //               else
    // 		    gfactor = 1.;
		  
    //               foutEps4 << x << " " << y << " " << 0.1973269718*gfactor*abs(lat->cells[pos]->getEpsilon()) << endl; 
    // 		  // abs just to get rid of negative 10^(-17) numbers at edge
    // 		}
    // 	      foutEps4 << endl;
    // 	    }
    // 	  foutEps4.close();
    // 	}

      if(it==1 || it==floor(it1) || it==floor(it2) || it==itmax)
	{	  
	  eccentricity(lat,group,param,it,0.1,0);
	  eccentricity(lat,group,param,it,1.,0);
	  eccentricity(lat,group,param,it,10.,0);
	}

 //      // output for hydro
//       if(it==itmax)
// 	{
// 	  int hx = param->getSizeOutput();
// 	  int hy = hx;
// 	  int heta = param->getEtaSizeOutput();
// 	  double hL = param->getLOutput();
// 	  double deta = param->getDetaOutput();

// 	  if(hL>L)
// 	    cout << "WARNING: hydro grid length larger than the computed one." << endl;
	  
// 	  int xpos, ypos, xposUp, yposUp;
// 	  double fracx, fracy, x1, x2;
// 	  double xlow, xhigh, ylow, yhigh;
// 	  int pos4;
// 	  double result, resultT00, resultT0x, resultT0y, resultT0eta;
// 	  double ha;
// 	  ha = hL/static_cast<double>(hx);
	  
	  
// 	  stringstream streH_name;
// 	  streH_name << "epsilonHydro" << param->getMPIRank() << ".dat";
// 	  string eH_name;
// 	  eH_name = streH_name.str();

// 	  ofstream foutEps2(eH_name.c_str(),ios::out); 
// 	  //	  ofstream foutEps3("epsilonPlot.dat",ios::out); 
// 	  //# <npart>= 355.155 etamax= 32 xmax= 131 ymax= 131 deta= 0.5 dx= 0.2 dy= 0.2
// 	  foutEps2 << "# dummy " << 1 << " etamax= " << heta
// 		   << " xmax= " << hx << " ymax= " << hy << " deta= " << deta 
// 		   << " dx= " << ha << " dy= " << ha << endl; 
	  
// 	  for(int ieta=0; ieta<heta; ieta++) // loop over all positions
// 	    {
// 	      for(int ix=0; ix<hx; ix++) // loop over all positions
// 		{
// 		  for(int iy=0; iy<hy; iy++)
// 		    {
// 		      x = -hL/2.+ha*ix;
// 		      y = -hL/2.+ha*iy;
		      
// 		      xpos = static_cast<int>(floor((x+L/2.)/a+0.0000000001));
// 		      ypos = static_cast<int>(floor((y+L/2.)/a+0.0000000001));
		      
// 		      if(xpos<N-1)
// 			xposUp = xpos+1;
// 		      else
// 			xposUp = xpos;
		      
// 		      if(ypos<N-1)
// 			yposUp = ypos+1;
// 		      else
// 			yposUp = ypos;
		      
// 		      xlow = -L/2.+a*xpos;
// 		      ylow = -L/2.+a*ypos;
		      
// 		      xhigh = -L/2.+a*xposUp;
// 		      yhigh = -L/2.+a*yposUp;
		      
// 		      fracx = (x-xlow)/ha;
		      
// 		      pos1 = xpos*N+ypos;
// 		      pos2 = xposUp*N+ypos;
// 		      pos3 = xpos*N+yposUp;
// 		      pos4 = xposUp*N+yposUp;
		      
// 		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
// 			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getEpsilon());
// 		      else
// 			x1 = 0.;
		      
// 		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
// 			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getEpsilon());
// 		      else
// 			x2 = 0.;
		      
// 		      fracy = (y-ylow)/ha;
		      
// 		      result = (1.-fracy)*x1+fracy*x2;
		

// 		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
// 			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getTtautau());
// 		      else
// 			x1 = 0.;
		      
// 		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
// 			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getTtautau());
// 		      else
// 			x2 = 0.;
		      
// 		      resultT00 = (1.-fracy)*x1+fracy*x2;
		
// 		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
// 			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getTtaux());
// 		      else
// 			x1 = 0.;
		      
// 		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
// 			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getTtaux());
// 		      else
// 			x2 = 0.;
		      
// 		      resultT0x = (1.-fracy)*x1+fracy*x2;
		
// 		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
// 			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getTtauy());
// 		      else
// 			x1 = 0.;
		      
// 		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
// 			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getTtauy());
// 		      else
// 			x2 = 0.;
		      
// 		      resultT0y = (1.-fracy)*x1+fracy*x2;
		
// 		      if(pos1>0 && pos1<(N)*(N) && pos2>0 && pos2<(N)*(N))
// 			x1 = (1-fracx)*abs(lat->cells[pos1]->getEpsilon())+fracx*abs(lat->cells[pos2]->getTtaueta());
// 		      else
// 			x1 = 0.;
		      
// 		      if(pos3>0 && pos3<N*N && pos4>0 && pos4<N*N)
// 			x2 = (1-fracx)*abs(lat->cells[pos3]->getEpsilon())+fracx*abs(lat->cells[pos4]->getTtaueta());
// 		      else
// 			x2 = 0.;
		      
// 		      resultT0eta = (1.-fracy)*x1+fracy*x2;
		
      
// 		      if(param->getRunningCoupling())
// 			{
// 			  if(pos>0 && pos<(N-1)*N+N-1)
// 			    {
// 			      g2mu2A = lat->cells[pos]->getg2mu2A();
// 			    }
// 			  else 
// 			    g2mu2A = 0;
			  
// 			  if(pos>0 && pos<(N-1)*N+N-1)
// 			    {
// 			      g2mu2B = lat->cells[pos]->getg2mu2B();
// 			    }
// 			  else
// 			    g2mu2B = 0;
			  
// 			  if(param->getRunWithQs()==2)
// 			    {
// 			      if(g2mu2A > g2mu2B)
// 				Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			      else
// 				Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			    }
// 			  else if(param->getRunWithQs()==0)
// 			    {
// 			      if(g2mu2A < g2mu2B)
// 				Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			      else
// 				Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			    }
// 			  else if(param->getRunWithQs()==1)
// 			    {
// 			      Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			    }
			  
// 			  if( param->getRunWithLocalQs() == 1 )
// 			    {
// 			      // 3 flavors
// 			      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
// 			      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			      
// 			      gfactor = g*g/(4.*PI*alphas);
// 			      //     cout << "Qs=" << Qs << endl;
// 			      //cout << alphas << endl;
// 			      // run with the local (in transverse plane) coupling
// 			    }
// 			  else
// 			    {
// 			      if ( param->getRunWithQs() == 0 )
// 				alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
// 			      else if ( param->getRunWithQs() == 1 )
// 				alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
// 			      else if ( param->getRunWithQs() == 2 )
// 				alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
			      
// 			      gfactor = g*g/(4.*PI*alphas);
// 			    }
// 			}
// 		      else
// 			gfactor = 1.;
		      		      
// 		      if( param->getWriteOutputs() == 1)
// 			foutEps2 << -heta/2.*deta+deta*ieta << " " << x << " " << y << " " << abs(0.1973269718*result*gfactor)
// 				 << " " << abs(0.1973269718*resultT00*gfactor) << " " << abs(0.1973269718*resultT0x*gfactor)
// 				 << " " << abs(0.1973269718*resultT0y*gfactor) 
// 				 << " " << abs(0.1973269718*resultT0eta*gfactor) << endl; 
		      
// // 		      if(ieta==0)
// // 			foutEps3 << -heta/2.*deta+deta*ieta << " " << x << " " << y << " " << abs(0.1973269718*result*gfactor) << endl; 
// 		    }
// 		  //		  Fouteps3 << endl;
// 		}
// 	      foutEps2 << endl;
// 	    }
// 	  foutEps2.close();
// 	  //	  foutEps3.close();
// 	}
//       // done output for hydro
      
      //      cout << "avgEps=" << avgEps << endl;

      //      foutEpsA << (it)*dtau*a << " " << 0.1973269718*avgEps*pow(a,2.) << endl;

      //" " << 0.1973269718*avgEpsMag*pow(a,2.) << " " << 0.1973269718*avgEpsEl*pow(a,2.) << " " << 0.1973269718*avgEpsBt*pow(a,2.) << " " << 0.1973269718*avgEpsEt*pow(a,2.) << endl;
      
      // ofstream foutdEdy(dEdy_name.c_str(),ios::app); 
      // foutdEdy << (it)*dtau*a << " " << 0.1973269718*(it)*dtau*a*a*a*avgEps << endl;
      // foutdEdy.close();

      // if(it==itmax)
      // 	{
      // 	  ofstream foutNE(NpartdEdy_name.c_str(),ios::app); 
      // 	  foutNE << param->getNpart() << " " <<  0.1973269718*(it)*dtau*a*a*a*avgEps << endl;
      // 	  foutNE.close();
      // 	}

      if(it==itmax)
	{	  
	  checkGaussLaw(lat, group, param, dtau, (it)*dtau);
	}
     
      // if( it % 10 == 0 && param->getWriteEvolution() )
      // 	{
      // 	  Tmunu(lat,group,param,it);
      // 	  anisotropy(lat,group,param,it);
      // 	}
      
      int success=1;
      if(it==1 || it==itmax)
	{	
	  success = multiplicity(lat,group,param,it);
	  //	  correlationsColor(lat,group,param,it);
	  //twoPointFunctionInK(param, lat, 0);
	}
      
      if (success==0)
	break;
    }
  //  foutEpsA.close();
  // foutEB.close();
}


void Evolution::Tmunu(Lattice *lat, Group *group, Parameters *param, int it)
{
  ofstream Tfout("Tmunu",ios::app); 

  double averageTtautau=0.;
  double averageTtaueta=0.;
  double averageTxx=0.;


  int N = param->getSize();
  int Nc = param->getNc(); 
  int pos, posNew, posX, posY, posmX, posmY, posmXmY, posXY, posmXpY, pospXmY, pos2X, pos2Y, posX2Y, pos2XY;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double Pi;
  Pi = param->getPi();
  double Ttautau=0.;
  double g = param->getg();
  double dtau = param->getdtau();
  double alphas = param->getalphas();
  double gfactor = g*g/(4.*3.141592654*alphas);
  double maxtime =  param->getMaxtime();
  Matrix one(Nc,1.);
 
  Matrix Ux(Nc);
  Matrix Uy(Nc);
  Matrix UxmX(Nc);
  Matrix UymY(Nc);
  Matrix UDx(Nc);
  Matrix UDy(Nc);
  Matrix UDxmX(Nc);
  Matrix UDymY(Nc);
  Matrix UDxmXpY(Nc);
  Matrix UDxpXpY(Nc);
  Matrix UxpX(Nc);
  Matrix UxpY(Nc);
  Matrix UDxpY(Nc);
  Matrix UxpXpY(Nc);
  Matrix UDypXmY(Nc);
  Matrix UypY(Nc);
  Matrix UypX(Nc);
  Matrix UDypX(Nc);
  Matrix UypXpY(Nc);
  Matrix UDypXpY(Nc);
  Matrix UymX(Nc);
  Matrix UxmXpY(Nc);
  Matrix UxmY(Nc);
  Matrix UDxmY(Nc);
  Matrix UypXmY(Nc);
  Matrix UDyp2X(Nc);
  Matrix Uyp2X(Nc);
  Matrix UDxpX(Nc);
  Matrix Uxp2Y(Nc);
  Matrix UDxp2Y(Nc);
  Matrix UDypY(Nc);
  Matrix UDymX(Nc);
  Matrix Uplaq(Nc), UplaqD(Nc), Uplaq1(Nc), Uplaq1D(Nc),Uplaq2(Nc),Uplaq3(Nc),Uplaq4(Nc);
  Matrix E1(Nc);
  Matrix E2(Nc);
  Matrix E1m(Nc);
  Matrix E2m(Nc);
  Matrix E1p(Nc);
  Matrix E2p(Nc);
  Matrix pi(Nc);
  Matrix piX(Nc);
  Matrix piY(Nc);
  Matrix piXY(Nc);
  Matrix phi(Nc);
  Matrix phiX(Nc);
  Matrix phiY(Nc);
  Matrix phiXY(Nc);
  Matrix phimX(Nc);
  Matrix phimY(Nc);
  Matrix phi2XY(Nc);
  Matrix phiX2Y(Nc);
  Matrix phi2X(Nc);
  Matrix phi2Y(Nc);
  Matrix phimXpY(Nc);
  Matrix phipXmY(Nc);
  Matrix phiTildeX(Nc);
  Matrix phiTildeY(Nc);
  Matrix phiTildeXY1(Nc);
  Matrix phiTildeXY2(Nc);
 

  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));

  //set plaquette in every cell
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;

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
	  lat->cells[pos]->setUplaq(Uplaq);
	  if(i==N-1)
	    lat->cells[pos]->setUplaq(one);
	}
    }
  
  //  stringstream strTmunu_name;
  //  strTmunu_name << "Tmunu" << param->getMPIRank();

  
  // just to make the result finite - doesn't make a difference for the components at tau=0
  // if (it==0)
  //   it = 1;


  // if(it==1)
  //   strTmunu_name << "Initial";
  // else if (it==itmax)
  //   strTmunu_name << "Final";
    
  // strTmunu_name << ".dat";

  // string Tmunu_name;
  // Tmunu_name = strTmunu_name.str();


  // ofstream foutTtt(Tmunu_name.c_str(),ios::out); 

  // T^\tau\tau, Txx, Tyy, Tetaeta:
  // electric part:
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  if(i<N-1)
	    posX = (i+1)*N+j;
	  else
	    posX = j;
	  if(j<N-1)
	    posY = i*N+(j+1);
	  else
	    posY = i*N;
	  
	  if(i<N-1 && j<N-1)
	    posXY =  (i+1)*N+j+1;
	  else if(i<N-1)
	    posXY =  (i+1)*N;
	  else if(j<N-1)
	    posXY = j+1;
	  else
	    posXY = 0;
	  
	  E1 = lat->cells[pos]->getE1();
	  E2 = lat->cells[pos]->getE2();
	  E1p = lat->cells[posY]->getE1(); 
	  E2p = lat->cells[posX]->getE2(); //shift y value in x direction
	
	  pi = lat->cells[pos]->getpi();
	  piX = lat->cells[posX]->getpi();
	  piY = lat->cells[posY]->getpi();
	  piXY = lat->cells[posXY]->getpi();

	  lat->cells[pos]->setTtautau((g*g/(it*dtau)/(it*dtau)*real((E1*E1).trace()+(E1p*E1p).trace()+(E2*E2).trace()+(E2p*E2p).trace())/2. //trans.
				       +(((pi*pi).trace()).real()+((piX*piX).trace()).real()
				        +((piY*piY).trace()).real()+((piXY*piXY).trace()).real())/4.
				       )); //long.
	  lat->cells[pos]->setTxx((g*g/(it*dtau)/(it*dtau)*real(-1.*(E1*E1).trace()-(E1p*E1p).trace()+(E2*E2).trace()+(E2p*E2p).trace())/2. //trans.
				       +(((pi*pi).trace()).real()+((piX*piX).trace()).real()
					 +((piY*piY).trace()).real()+((piXY*piXY).trace()).real())/4.)); //long.
	  lat->cells[pos]->setTyy((g*g/(it*dtau)/(it*dtau)*real((E1*E1).trace()+(E1p*E1p).trace()-(E2*E2).trace()-(E2p*E2p).trace())/2. //trans.
				       +(((pi*pi).trace()).real()+((piX*piX).trace()).real()
					 +((piY*piY).trace()).real()+((piXY*piXY).trace()).real())/4.)); //long.
	  lat->cells[pos]->setTetaeta(1./(it*dtau)/(it*dtau)*
				      ((g*g/(it*dtau)/(it*dtau)*real((E1*E1).trace()+(E1p*E1p).trace()+(E2*E2).trace()+(E2p*E2p).trace())/2. //trans.
				       -(((pi*pi).trace()).real()+((piX*piX).trace()).real()
					 +((piY*piY).trace()).real()+((piXY*piXY).trace()).real())/4.))); //long.
	}
    } 
  // magnetic part:
  stringstream strT_name;
  strT_name << "T" << param->getMPIRank() << ".dat";
  string T_name;
  T_name = strT_name.str();

  ofstream fout(T_name.c_str(),ios::out); 
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  if(i<N-1)
	    posX = (i+1)*N+j;
	  else
	    posX = j;
	  if(j<N-1)
	    posY = i*N+(j+1);
	  else
	    posY = i*N;
	  
	  if(i<N-1 && j<N-1)
	    posXY =  (i+1)*N+j+1;
	  else if(i<N-1)
	    posXY =  (i+1)*N;
	  else if(j<N-1)
	    posXY = j+1;
	  else
	    posXY = 0;
	  
	  Uplaq = lat->cells[pos]->getUplaq();
	  
	  phi = lat->cells[pos]->getphi();
	  phiX = lat->cells[posX]->getphi();
	  phiY = lat->cells[posY]->getphi();
	  phiXY = lat->cells[posXY]->getphi();
	
	  Ux = lat->cells[pos]->getUx();
	  Uy = lat->cells[pos]->getUy();
	  UDx = Ux;
	  UDx.conjg();
	  UDy = Uy;
	  UDy.conjg();
	  
	  phiTildeX = Ux*phiX*UDx;
	  phiTildeY = Uy*phiY*UDy;

	  //same at one up in he other direction
	  Ux = lat->cells[posY]->getUx();
	  Uy = lat->cells[posX]->getUy();
	  UDx = Ux;
	  UDx.conjg();
	  UDy = Uy;
	  UDy.conjg();
	  
	  phiTildeXY1 = Ux*phiXY*UDx;
	  phiTildeXY2 = Uy*phiXY*UDy;
	  
	  lat->cells[pos]->setTtautau(
				      lat->cells[pos]->getTtautau() 
				      + 2./pow(g,2.)*(static_cast<double>(Nc)-(Uplaq.trace()).real())
				      +
				      0.5/(it*dtau)/(it*dtau)*(real(((phi-phiTildeX)*(phi-phiTildeX)).trace())
							       +real(((phiY-phiTildeXY1)*(phiY-phiTildeXY1)).trace())
							       +real(((phi-phiTildeY)*(phi-phiTildeY)).trace())
				       			       +real(((phiX-phiTildeXY2)*(phiX-phiTildeXY2)).trace()))
				      );

	  lat->cells[pos]->setTxx(
				  lat->cells[pos]->getTxx()
				  + 2./pow(g,2.)*(static_cast<double>(Nc)-(Uplaq.trace()).real())
				  +
				  0.5/(it*dtau)/(it*dtau)*(real(((phi-phiTildeX)*(phi-phiTildeX)).trace())
							   +real(((phiY-phiTildeXY1)*(phiY-phiTildeXY1)).trace())
							   -real(((phi-phiTildeY)*(phi-phiTildeY)).trace())
							   -real(((phiX-phiTildeXY2)*(phiX-phiTildeXY2)).trace()))
				  );

	  lat->cells[pos]->setTyy(
				  lat->cells[pos]->getTyy()
				  + 2./pow(g,2.)*(static_cast<double>(Nc)-(Uplaq.trace()).real())
				  +
				  0.5/(it*dtau)/(it*dtau)*(-real(((phi-phiTildeX)*(phi-phiTildeX)).trace())
							   -real(((phiY-phiTildeXY1)*(phiY-phiTildeXY1)).trace())
							   +real(((phi-phiTildeY)*(phi-phiTildeY)).trace())
							   +real(((phiX-phiTildeXY2)*(phiX-phiTildeXY2)).trace()))
				  );

	  lat->cells[pos]->setTetaeta(
				      lat->cells[pos]->getTetaeta()
				      + 1./(it*dtau)/(it*dtau)*(
								-2./pow(g,2.)*(Nc-(Uplaq.trace()).real())
								+
								0.5/(it*dtau)/(it*dtau)*(+real(((phi-phiTildeX)*(phi-phiTildeX)).trace())
											 +real(((phiY-phiTildeXY1)*(phiY-phiTildeXY1)).trace())
											 +real(((phi-phiTildeY)*(phi-phiTildeY)).trace())
											 +real(((phiX-phiTildeXY2)*(phiX-phiTildeXY2)).trace()))
								)
				      );

	  // clean up numerical noise outside the interaction region
	  if (lat->cells[pos]->getg2mu2A()<1e-12 || lat->cells[pos]->getg2mu2B()<1e-12)
	    lat->cells[pos]->setEpsilon(0.);
	  else
	    lat->cells[pos]->setEpsilon(lat->cells[pos]->getTtautau()*1/pow(a,4.));

	  
	  lat->cells[pos]->setTtautau(lat->cells[pos]->getTtautau()*1/pow(a,4.));
	  lat->cells[pos]->setTxx(lat->cells[pos]->getTxx()*1/pow(a,4.));
	  lat->cells[pos]->setTyy(lat->cells[pos]->getTyy()*1/pow(a,4.));
	  lat->cells[pos]->setTetaeta(lat->cells[pos]->getTetaeta()*1/pow(a,6.));

	  //	  foutTtt << i << " " << j << " " << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTxx();
	  //foutTtt << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTetaeta() << "\n";
	}
      //      foutTtt << endl;
    }
  
  fout.close();
  //  foutTtt.close();

  // T^\tau x, T^\tau y
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  if(i<N-1)
	    posX = (i+1)*N+j;
	  else
	    posX = j;
	  if(j<N-1)
	    posY = i*N+(j+1);
	  else
	    posY = i*N;
	  
	  if(i<N-1 && j<N-1)
	    posXY =  (i+1)*N+j+1;
	  else if(i<N-1)
	    posXY =  (i+1)*N;
	  else if(j<N-1)
	    posXY = j+1;
	  else
	    posXY = 0;

	  if(i>0)
	    posmX = (i-1)*N+j;
	  else
	    posmX = (N-1)*N+j;

	  if(j>0)
	    posmY = i*N+(j-1);
	  else
	    posmY = i*N+N-1;

	  if(i>0 && j<N-1)
	    posmXpY =  (i-1)*N+j+1;
	  else if(i>0)
	    posmXpY =  (i-1)*N;
	  else if(j<N-1)
	    posmXpY = (N-1)*N+j+1;
	  else
	    posmXpY = (N-1)*N;
	  
	  if(j>0 && i<N-1)
	    pospXmY =  (i+1)*N+j-1;
	  else if(j>0)
	    pospXmY =  j-1;
	  else if(i<N-1)
	    pospXmY =  (i+1)*N+N-1;
	  else
	    pospXmY = N-1;

	  if(i<N-2)
	    pos2X = (i+2)*N+j;
	  else
	    pos2X = (i+2-N)*N+j;
	  if(j<N-2)
	    pos2Y = i*N+(j+2);
	  else
	    pos2Y = i*N+(j+2-N);

	  if(i<N-2 && j<N-1)
	    pos2XY =  (i+2)*N+j+1;
	  else if(i<N-2)
	    pos2XY =  (i+2)*N;
	  else if(j<N-1)
	    pos2XY = (i+2-N)*N+j+1;
	  else
	    pos2XY = (i+2-N)*N;

	  if(i<N-1 && j<N-2)
	    posX2Y =  (i+1)*N+j+2;
	  else if(i<N-1)
	    posX2Y =  (i+1)*N+j+2-N;
	  else if(j<N-2)
	    posX2Y = j+2;
	  else
	    posX2Y = j+2-N;

	  E1 = lat->cells[pos]->getE1();
	  E2 = lat->cells[pos]->getE2();
	  E1p = lat->cells[posY]->getE1(); //shift x value in y direction
	  E2p = lat->cells[posX]->getE2(); // shift y value in x direction

	  Uplaq = lat->cells[pos]->getUplaq();
	  Uplaq1 = lat->cells[posmX]->getUplaq();
	  Uplaq2 = lat->cells[posmY]->getUplaq();
	  UplaqD = Uplaq; 
	  UplaqD.conjg();
	  Uplaq1D = Uplaq1; 
	  Uplaq1D.conjg();
	  
	  pi = lat->cells[pos]->getpi();
	  piX = lat->cells[posX]->getpi();
	  piY = lat->cells[posY]->getpi();
	  piXY = lat->cells[posXY]->getpi();

	  phi = lat->cells[pos]->getphi();
	  phimX = lat->cells[posmX]->getphi();
	  phiX = lat->cells[posX]->getphi();
	  phimY = lat->cells[posmY]->getphi();
	  phiY = lat->cells[posY]->getphi();
	  phiXY = lat->cells[posXY]->getphi();
	  phimXpY = lat->cells[posmXpY]->getphi();
	  phipXmY = lat->cells[pospXmY]->getphi();
	  phi2X = lat->cells[pos2X]->getphi();
	  phi2XY = lat->cells[pos2XY]->getphi();
	  phi2Y = lat->cells[pos2Y]->getphi();
	  phiX2Y = lat->cells[posX2Y]->getphi();

	  Ux = lat->cells[pos]->getUx();
	  UDx = Ux;
	  UDx.conjg();
	
	  UxmX = lat->cells[posmX]->getUx();
	  UDxmX = lat->cells[posmX]->getUx();
	  UDxmX.conjg();
	  UxmXpY = lat->cells[posmXpY]->getUx();
	  UDxmXpY = lat->cells[posmXpY]->getUx();
	  UDxmXpY.conjg();
	  
	  UxpX = lat->cells[posX]->getUx();
	  UxpY = lat->cells[posY]->getUx();
	  UDxpX = UxpX;
	  UDxpX.conjg();
	  UDxpY = UxpY;
	  UDxpY.conjg();

	  UxpXpY = lat->cells[posXY]->getUx();
	  UDxpXpY = lat->cells[posXY]->getUx();
	  UDxpXpY.conjg();
	  UxmXpY = lat->cells[posmXpY]->getUx();
	  UDxmXpY = UxmXpY;
	  UDxmXpY.conjg();

	  Uy = lat->cells[pos]->getUy();
	  UDy = Uy;
	  UDy.conjg();
	  
	  UymY = lat->cells[posmY]->getUy();
	  UDymY = lat->cells[posmY]->getUy();
	  UDymY.conjg();
	  UypXmY = lat->cells[pospXmY]->getUy();
	  UDypXmY = lat->cells[pospXmY]->getUy();
	  UDypXmY.conjg();
	  
	  UypY = lat->cells[posY]->getUy();
	  UypX = lat->cells[posX]->getUy();
	  UDypX = UypX;
	  UDypX.conjg();

	  UDypY = lat->cells[posY]->getUy();
	  UDypY.conjg();

	  UDxpX = lat->cells[posX]->getUx();
	  UDxpX.conjg();

	  Uyp2X = lat->cells[pos2X]->getUy();
	  UDyp2X = Uyp2X;
	  UDyp2X.conjg();
	  Uxp2Y = lat->cells[pos2Y]->getUx();
	  UDxp2Y = Uxp2Y;
	  UDxp2Y.conjg();
	  
	  UypXpY = lat->cells[posXY]->getUy();
	  UDypXpY = lat->cells[posXY]->getUy();
	  UDypXpY.conjg();
	  UymX = lat->cells[posmX]->getUy();
	  UDymX = UymX;
	  UDymX.conjg();
	  UxmY = lat->cells[posmY]->getUx();
	  UDxmY = UxmY;
	  UDxmY.conjg();
	  UypXmY = lat->cells[pospXmY]->getUy();
	  


// 	  lat->cells[pos]->setTtaux( 
// 				    +2./(it*dtau)/4. * (E2*(UxpY*UDypX*UDx*Uy+UDy*UDxmX*UymX*UxmXpY)
// 							+E2p*(UxpXpY*UDyp2X*UDxpX*UypX+UDypX*UDx*Uy*UxpY)).trace().imag()
// 				    -2./8./(it*dtau) * ( 
// 					     		pi * (Ux*phiX*UDx - UDxmX*phimX*UxmX)
// 							+ piY * (UxpY*phiXY*UDxpY - UDxmXpY*phimXpY*UxmXpY)
// 							+ piX * (UxpX*phi2X*UDxpX - UDx*phi*Ux)
// 							+ piXY * (UxpXpY*phi2XY*UDxpXpY - UDxpY*phiY*UxpY)
// 							 ).trace().real()
// 				     );



// 	  lat->cells[pos]->setTtauy( 
// 				    +2./(it*dtau)/4. * (E1*(UypX*UDxpY*UDy*Ux+UDx*UDymY*UxmY*UypXmY)
// 							+E1p*(UypXpY*UDxp2Y*UDypY*UxpY+UDxpY*UDy*Ux*UypX)).trace().imag()
// 				    -2./8./(it*dtau) * ( 
// 							pi * (Uy*phiY*UDy - UDymY*phimY*UymY)
// 							+ piX * (UypX*phiXY*UDypX - UDypXmY*phipXmY*UypXmY)
// 							+ piY * (UypY*phi2Y*UDypY - UDy*phi*Uy)
// 							+ piXY * (UypXpY*phiX2Y*UDypXpY - UDypX*phiX*UypX)
// 							 ).trace().real()
// 				     );


	  // note that the minus sign of the first terms in T^\taux and T^\tauy comes from the direction of the plaquettes - I am using +F^{yx} instead
	  // of -F^{xy} if you like.
	  lat->cells[pos]->setTtaux( 
				    +2./(it*dtau)/8. * (E2*(Uy*UxpY*UDypX*UDx - Ux*UypX*UDxpY*UDy
							    - (Uy*UxpY*UDypX*UDx-Ux*UypX*UDxpY*UDy).trace()/static_cast<double>(Nc)*one
							    +UDxmX*UymX*UxmXpY*UDy - Uy*UDxmXpY*UDymX*UxmX
							    -(UDxmX*UymX*UxmXpY*UDy - Uy*UDxmXpY*UDymX*UxmX).trace()/static_cast<double>(Nc)*one)
							+E2p*(UypX*UxpXpY*UDyp2X*UDxpX-UxpX*Uyp2X*UDxpXpY*UDypX
							      -(UypX*UxpXpY*UDyp2X*UDxpX-UxpX*Uyp2X*UDxpXpY*UDypX).trace()/static_cast<double>(Nc)*one
							      +UDx*Uy*UxpY*UDypX-UypX*UDxpY*UDy*Ux
							      -(UDx*Uy*UxpY*UDypX-UypX*UDxpY*UDy*Ux).trace()/static_cast<double>(Nc)*one)
							).trace().imag()
				    -2./8./(it*dtau) * ( 
					     		pi * (Ux*phiX*UDx - UDxmX*phimX*UxmX)
							+ piY * (UxpY*phiXY*UDxpY - UDxmXpY*phimXpY*UxmXpY)
							+ piX * (UxpX*phi2X*UDxpX - UDx*phi*Ux)
							+ piXY * (UxpXpY*phi2XY*UDxpXpY - UDxpY*phiY*UxpY)
							 ).trace().real()
				     );



	  lat->cells[pos]->setTtauy( 
				    +2./(it*dtau)/8. * (E1*(Ux*UypX*UDxpY*UDy-Uy*UxpY*UDypX*UDx
							    -(Ux*UypX*UDxpY*UDy-Uy*UxpY*UDypX*UDx).trace()/static_cast<double>(Nc)*one
							    +UDymY*UxmY*UypXmY*UDx-Ux*UDypXmY*UDxmY*UymY
							    -(UDymY*UxmY*UypXmY*UDx-Ux*UDypXmY*UDxmY*UymY).trace()/static_cast<double>(Nc)*one
							    )
							+E1p*(UxpY*UypXpY*UDxp2Y*UDypY-UypY*Uxp2Y*UDypXpY*UDxpY
							      -(UxpY*UypXpY*UDxp2Y*UDypY-UypY*Uxp2Y*UDypXpY*UDxpY).trace()/static_cast<double>(Nc)*one
							      +UDy*Ux*UypX*UDxpY-UxpY*UDypX*UDx*Uy
							      -(UDy*Ux*UypX*UDxpY-UxpY*UDypX*UDx*Uy).trace()/static_cast<double>(Nc)*one
							      )
							).trace().imag()
				    -2./8./(it*dtau) * ( 
							pi * (Uy*phiY*UDy - UDymY*phimY*UymY)
							+ piX * (UypX*phiXY*UDypX - UDypXmY*phipXmY*UypXmY)
							+ piY * (UypY*phi2Y*UDypY - UDy*phi*Uy)
							+ piXY * (UypXpY*phiX2Y*UDypXpY - UDypX*phiX*UypX)
							 ).trace().real()
				     );





	  lat->cells[pos]->setTtaueta( 
				      g/(it*dtau)/(it*dtau)/(it*dtau)* 
				      ( 
				       E1 * (Ux*phiX*UDx - phi)
				       + E1p * (UxpY*phiXY*UDxpY - phiY)
				       + E2 * (Uy*phiY*UDy - phi)
				       + E2p * (UypX*phiXY*UDypX - phiX)
					).trace().real()
				       );
	  

	  // T^xy
	  
	  lat->cells[pos]->setTxy(2./(it*dtau)/(it*dtau) * (
							    -1./4.*g*g*(E1+Uy*E1p*UDy)*(E2+Ux*E2p*UDx)
							    //-1./4.*g*g*(E1+E1p)*(E2+E2p)
							    +1./4.*( 
								    (Ux*phiX*UDx-phi)*(Uy*phiY*UDy-phi) 
								     + Uy*(UxpY*phiXY*UDxpY - phiY)*UDy*(Uy*phiY*UDy-phi)
								     + (Ux*phiX*UDx - phi)*Ux*(UypX*phiXY*UDypX - phiX)*UDx
								     + Uy*(UxpY*phiXY*UDxpY - phiY)*UDy*Ux*(UypX*phiXY*UDypX - phiX)*UDx
								     )
			       				    ).trace().real()
				  ); 

	  lat->cells[pos]->setTxeta(-2./(it*dtau)/(it*dtau) * (1./4.*g*(E1*(pi+Ux*piX*UDx)+E1p*(piY+UxpY*piXY*UDxpY)).trace().real()
							       +1./8./g*( (Ux*UypX*UDxpY*UDy-Uy*UxpY*UDypX*UDx
									 -(Ux*UypX*UDxpY*UDy-Uy*UxpY*UDypX*UDx).trace()/static_cast<double>(Nc)*one
									 + Uy*UDxmXpY*UDymX*UxmX-UDxmX*UymX*UxmXpY*UDy
									 -(Uy*UDxmXpY*UDymX*UxmX-UDxmX*UymX*UxmXpY*UDy).trace()/static_cast<double>(Nc)*one)
									*(Uy*phiY*UDy-phi)
									  +
									(UypX*UDxpY*UDy*Ux-UDx*Uy*UxpY*UDypX
									 -(UypX*UDxpY*UDy*Ux-UDx*Uy*UxpY*UDypX).trace()/static_cast<double>(Nc)*one
									 +UxpX*Uyp2X*UDxpXpY*UDypX-UypX*UxpXpY*UDyp2X*UDxpX
									 -(UxpX*Uyp2X*UDxpXpY*UDypX-UypX*UxpXpY*UDyp2X*UDxpX).trace()/static_cast<double>(Nc)*one)
									*(UypX*phiXY*UDypX-phiX)
									).trace().imag()
							       )
				    );

	  lat->cells[pos]->setTyeta(-2./(it*dtau)/(it*dtau) * (1./4.*g*(E2*(pi+Uy*piY*UDy)+E2p*(piX+UypX*piXY*UDypX)).trace().real()
							       +1./8./g*( (Uy*UxpY*UDypX*UDx-Ux*UypX*UDxpY*UDy
									   -(Uy*UxpY*UDypX*UDx-Ux*UypX*UDxpY*UDy).trace()/static_cast<double>(Nc)*one
									   +Ux*UDypXmY*UDxmY*UymY-UDymY*UxmY*UypXmY*UDx
									   -(Ux*UDypXmY*UDxmY*UymY-UDymY*UxmY*UypXmY*UDx).trace()/static_cast<double>(Nc)*one)
									  *(Ux*phiX*UDx-phi)
									  +
									  (UxpY*UDypX*UDx*Uy-UDy*Ux*UypX*UDxpY
									   -(UxpY*UDypX*UDx*Uy-UDy*Ux*UypX*UDxpY).trace()/static_cast<double>(Nc)*one
									   +UypY*Uxp2Y*UDypXpY*UDxpY-UxpY*UypXpY*UDxp2Y*UDypY
									   -(UypY*Uxp2Y*UDypXpY*UDxpY-UxpY*UypXpY*UDxp2Y*UDypY).trace()/static_cast<double>(Nc)*one)
									  *(UxpY*phiXY*UDxpY-phiY)
									  ).trace().imag()
							       )
				    );
	  
	  

// 	  if(i==N/2 &&j==N/2)
// 	    {
// 	      cout << "normal E product:" << (E1+E1p)*(E2+E2p) << endl << endl;
// 	      cout << "U E product:" << (E1*E2+E1*Ux*E2p+E1p*UDy*E2+E1p*UDy*Ux*E2p) << endl << endl;
// 	    }


	  lat->cells[pos]->setTtaux(lat->cells[pos]->getTtaux()*1/pow(a,4.));
	  lat->cells[pos]->setTtauy(lat->cells[pos]->getTtauy()*1/pow(a,4.));
	  lat->cells[pos]->setTtaueta(lat->cells[pos]->getTtaueta()*1/pow(a,5.));
	  lat->cells[pos]->setTxy(lat->cells[pos]->getTxy()*1/pow(a,4.));
	  lat->cells[pos]->setTxeta(lat->cells[pos]->getTxeta()*1/pow(a,5.));
	  lat->cells[pos]->setTyeta(lat->cells[pos]->getTyeta()*1/pow(a,5.));





// 	  if(abs(lat->cells[pos]->getTtaux())<1e-7)
	 // 	    lat->cells[pos]->setTtaux(0.);
// 	  if(abs(lat->cells[pos]->getTtauy())<5)
// 	    lat->cells[pos]->setTtauy(0.);
// 	  if(abs(lat->cells[pos]->getTxy())<5)
// 	    lat->cells[pos]->setTxy(0.);
	  
	  
	  averageTtautau += lat->cells[pos]->getTtautau()*lat->cells[pos]->getTtautau();
	  averageTtaueta += lat->cells[pos]->getTtaueta()*lat->cells[pos]->getTtaueta();
	  averageTxx += lat->cells[pos]->getTxx()*lat->cells[pos]->getTxx();
	  
	  // if((i==300 && j==300))
  	  //   {
  	  //     cout << j << endl << endl;
  	  //     cout << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy()  << " " << lat->cells[pos]->getTtaueta() << endl;
  	  //     cout << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTxeta() << endl;
  	  //     cout << lat->cells[pos]->getTtauy() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta() << endl;
  	  //     cout << lat->cells[pos]->getTtaueta() << " " << lat->cells[pos]->getTxeta() << " " << lat->cells[pos]->getTyeta() << " " << lat->cells[pos]->getTetaeta() << endl;
  	  //     cout << endl;

	  //     //	      Tfout << it*dtau*a << " " << lat->cells[pos]->getTtautau() << " " << lat->cells[pos]->getTtaux() << " " << lat->cells[pos]->getTtauy()  << " " << lat->cells[pos]->getTtaueta() << " " << lat->cells[pos]->getTxx() << " " << lat->cells[pos]->getTxy() << " " << lat->cells[pos]->getTxeta() << " " << lat->cells[pos]->getTyy() << " " << lat->cells[pos]->getTyeta() << " " << lat->cells[pos]->getTetaeta() << endl;
  	  //   }
	}
    }
  averageTtautau /= double(N);
  averageTtaueta /= double(N);
  averageTxx /= double(N);
  Tfout << it*dtau*a << " " << sqrt(averageTtautau) << " " << sqrt(averageTtaueta) << " " << sqrt(averageTxx) << endl;
  Tfout.close();
}

void Evolution::u(Lattice *lat, Group *group, Parameters *param, int it)
{
  MyEigen *myeigen;
  myeigen = new MyEigen();
  
  myeigen->flowVelocity4D(lat,group,param,it);
  //cout << "3D:" << endl;
  //eigen->flowVelocity(lat,group,param,it);

  delete myeigen;
}


void Evolution::anisotropy(Lattice *lat, Group *group, Parameters *param, int it)
{
  stringstream straniso_name;
  straniso_name << "anisotropy" << param->getMPIRank() << ".dat";
  string aniso_name;
  aniso_name = straniso_name.str();

  ofstream foutAniso(aniso_name.c_str(),ios::app); 
  int N = param->getSize();
  double L = param->getL();
  double a = L/N; // lattice spacing in fm

  double num=0., den=0.;
  int pos;
  for(int ix=0; ix<N; ix++) 
    {
      for(int iy=0; iy<N; iy++)
	{
	  pos = ix*N+iy;
	  if (lat->cells[pos]->getTtautau() > 10.
	      )
	    {
	      num += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
	      den += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
	    }
	}
    }
  
  foutAniso << it*a*param->getdtau() << " " << num/den << endl; 

  foutAniso.close();
}


void Evolution::eccentricity(Lattice *lat, Group *group, Parameters *param, int it, double cutoff, int doAniso)
{
  stringstream strecc_name;
  strecc_name << "eccentricities" << param->getMPIRank() << ".dat";
  string ecc_name;
  ecc_name = strecc_name.str();
 
  // cutoff on energy density is 'cutoff' times Lambda_QCD^4
  int N = param->getSize();
  int pos, posNew;
  double rA, phiA, x, y;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double eccentricity1, eccentricity2, eccentricity3, eccentricity4, eccentricity5, eccentricity6;
  double avcos, avsin, avcos1, avsin1, avcos3, avsin3, avrSq, avxSq, avySq, avr1, avr3, avcos4, avsin4, avr4, avcos5, avsin5, avr5, avcos6, avsin6, avr6;
  double Rbar;
  double Pi;
  Pi = param->getPi();
  double Psi1, Psi2, Psi3, Psi4, Psi5, Psi6;
  double maxEps = 0;
  double g = param->getg();

  double g2mu2A, g2mu2B, gfactor, alphas, Qs;
  double c = param->getc();
  double muZero = param->getMuZero();

  double weight;

  double area = 0.;
  
  avrSq=0.;
  avr3=0.;

  double avx = 0.;
  double avy = 0.;
  double toteps = 0.;
  int xshift;
  int yshift;
  double maxX = 0.;
  double maxY = 0.;

  double smallestX = 0.;
  double smallestY = 0.;
  
  for(int ix=0; ix<N; ix++) 
    {
      for(int iy=0; iy<N; iy++)
	{
	  pos = ix*N+iy;
	  maxEps = max(lat->cells[pos]->getEpsilon(),maxEps);
	}
    }

  // cout << "maxEps=" << maxEps << endl;

  // first shift to the center
  for(int ix=0; ix<N; ix++) 
    {
      x = -L/2.+a*ix;
      for(int iy=0; iy<N; iy++)
	{
	  y = -L/2.+a*iy;
	  pos = ix*N+iy;
   
	      if(param->getRunningCoupling())
		{
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2A = lat->cells[pos]->getg2mu2A();
		    }
		  else 
		    g2mu2A = 0;
		  
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2B = lat->cells[pos]->getg2mu2B();
		    }
		  else
		    g2mu2B = 0;
		  
		  if(param->getRunWithQs()==2)
		    {
		      if(g2mu2A > g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==0)
		    {
		      if(g2mu2A < g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==1)
		    {
		      Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  
		  if ( param->getRunWithLocalQs() == 1 )
		    {
		      // 3 flavors
		      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		      
		      gfactor = g*g/(4.*PI*alphas);
		      //     cout << "Qs=" << Qs << endl;
		      //cout << alphas << endl;
		      // run with the local (in transverse plane) coupling
		    }
		  else
		    {
		      if ( param->getRunWithQs() == 0 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 1 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 2 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		      
		      gfactor = g*g/(4.*PI*alphas);
		    }
		}
	      else
		gfactor = 1.;	      
	      
 	      if (lat->cells[pos]->getEpsilon()*gfactor < cutoff)  //this is 1/fm^4, so Lambda_QCD^{-4} (because \Lambda_QCD is roughly 1/fm)
		{
		  // if (lat->cells[pos]->getEpsilon()*gfactor >0.)
		  //   cout << lat->cells[pos]->getEpsilon()*gfactor  << endl;
		  weight = 0.;
		}
	      else
		{
		  weight = lat->cells[pos]->getEpsilon()*gfactor;
		  area += a*a;
		}
	      avx += x*weight;
	      avy += y*weight;
	      toteps += weight;    
	}
    }
      
  avx/=toteps;
  avy/=toteps;
  
  param->setArea(0.);
  if(it==1 && cutoff==10)
    param->setArea(area);

  //  cout << "avx=" << avx << endl;
  //cout << "avy=" << avy << endl;
  
  xshift = static_cast<int>(floor(avx/a+0.00000000001));
  yshift = static_cast<int>(floor(avy/a+0.00000000001));
  
  //  cout << "xshift=" << xshift << endl;
  //cout << "yshift=" << yshift << endl;
  
  //  cout << "a xshift=" << a*xshift << endl;
  //cout << "a yshift=" << a*yshift << endl;
  
  avcos1 = 0.;
  avsin1 = 0.;
  avcos = 0.;
  avsin = 0.;
  avcos3 = 0.;
  avsin3 = 0.;
  avcos4 = 0.;
  avsin4 = 0.;
  avcos5 = 0.;
  avsin5 = 0.;
  avcos6 = 0.;
  avsin6 = 0.;
  avr1=0.;
  avrSq=0.;
  avxSq=0.;
  avySq=0.;
  avr3=0.;
  avr4=0.;
  avr5=0.;
  avr6=0.;

   
  for(int ix=2; ix<N-2; ix++) 
	{
	  x = -L/2.+a*ix-avx;
	  for(int iy=2; iy<N-2; iy++)
	    {
	      pos = ix*N+iy;
	      y = -L/2.+a*iy-avy;
	      if (x>=0)
		{
		  phiA = atan(y/x);
		  if (x==0) 
		    {
		      if (y>=0) phiA=PI/2.;
		      else if (y<0) phiA=3.*PI/2.;
		    }
		}
	      else
		{
		  phiA = atan(y/x)+PI;
		}

 	      // check this
//  	      if(lat->cells[pos]->getEpsilon()>maxEps/100000.)
//  		cout << x << " " << y << " " << lat->cells[pos]->getEpsilon() << endl;

	      
	      if(param->getRunningCoupling())
		{
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2A = lat->cells[pos]->getg2mu2A();
		    }
		  else 
		    g2mu2A = 0;
		  
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2B = lat->cells[pos]->getg2mu2B();
		    }
		  else
		    g2mu2B = 0;
		  
		  if(param->getRunWithQs()==2)
		    {
		      if(g2mu2A > g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==0)
		    {
		      if(g2mu2A < g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==1)
		    {
		      Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  
		  if ( param->getRunWithLocalQs() == 1 )
		    {
		      // 3 flavors
		      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		      
		      gfactor = g*g/(4.*PI*alphas);
		      //     cout << "Qs=" << Qs << endl;
		      //cout << alphas << endl;
		      // run with the local (in transverse plane) coupling
		    }
		  else
		    {
		      if ( param->getRunWithQs() == 0 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 1 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 2 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		      
		      gfactor = g*g/(4.*PI*alphas);
		    }
		}
	      else
		gfactor = 1.;
	      
	      
 	      if (lat->cells[pos]->getEpsilon()*gfactor < cutoff)  //this is 1/fm^4, so Lambda_QCD^{-4}
		{
		  weight = 0.;
		}
	      else
		{
		  weight = lat->cells[pos]->getEpsilon()*gfactor;
		}
	      
	      rA = sqrt( x*x + y*y );
	      avr1 += rA*rA*rA*(weight);
	      avrSq += rA*rA*(weight); // compute average r^2
	      avr3 += rA*rA*rA*(weight);
	      avr4 += rA*rA*rA*rA*(weight);
	      avr5 += rA*rA*rA*rA*rA*(weight);
	      avr6 += rA*rA*rA*rA*rA*rA*(weight);
	      
	      avcos1 += rA*rA*rA*cos(phiA)*(weight);
	      avsin1 += rA*rA*rA*sin(phiA)*(weight);
	      avcos  += rA*rA*cos(2.*phiA)*(weight);
	      avsin  += rA*rA*sin(2.*phiA)*(weight);
	      avcos3 += rA*rA*rA*cos(3.*phiA)*(weight);
	      avsin3 += rA*rA*rA*sin(3.*phiA)*(weight);
	      avcos4 += rA*rA*rA*rA*cos(4.*phiA)*(weight);
	      avsin4 += rA*rA*rA*rA*sin(4.*phiA)*(weight);
	      avcos5 += rA*rA*rA*rA*rA*cos(5.*phiA)*(weight);
	      avsin5 += rA*rA*rA*rA*rA*sin(5.*phiA)*(weight);
	      avcos6 += rA*rA*rA*rA*rA*rA*cos(6.*phiA)*(weight);
	      avsin6 += rA*rA*rA*rA*rA*rA*sin(6.*phiA)*(weight);
	      
	      if(weight > cutoff && iy == N/2 + yshift )
		{
		  maxX = x;
		}
	      if(weight > cutoff && ix == N/2 + xshift)
		{
		  maxY = y;
		}

	      if(weight < cutoff && iy == N/2 + yshift && ix > N/2 + xshift && smallestX==0)
		{
		  smallestX = x;
		}
	      if(weight < cutoff && ix == N/2 + xshift && iy > N/2 + yshift  && smallestY==0)
		{
		  smallestY = y;
		}
	    }
	}
      
      // compute and print eccentricity and angles:
      Psi1 = (atan(avsin1/avcos1)+PI)/1.;
      Psi2 = (atan(avsin/avcos)+PI)/2.;
      Psi3 = (atan(avsin3/avcos3)+PI)/3.;
      Psi4 = (atan(avsin4/avcos4)+PI)/4.;
      Psi5 = (atan(avsin5/avcos5)+PI)/5.;
      Psi6 = (atan(avsin6/avcos6)+PI)/6.;
      eccentricity1 = sqrt(avcos1*avcos1+avsin1*avsin1)/avr1;
      eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
      eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avr3;
      eccentricity4 = sqrt(avcos4*avcos4+avsin4*avsin4)/avr4;
      eccentricity5 = sqrt(avcos5*avcos5+avsin5*avsin5)/avr5;
      eccentricity6 = sqrt(avcos6*avcos6+avsin6*avsin6)/avr6;
   
      
      //      cout << "ecc1=" << eccentricity1 << endl;
      //cout << "ecc2=" << eccentricity2 << endl;
      //cout << "ecc3=" << eccentricity3 << endl;
      //cout << "ecc4=" << eccentricity4 << endl;
      //cout << "ecc5=" << eccentricity5 << endl;
      //cout << "ecc6=" << eccentricity6 << endl;
   
      double avx2 = avx;
      double avy2 = avy;
      avx=0.;
      avy=0.;
      toteps=0.;

      for(int ix=0; ix<N; ix++) 
	{
	  x = -L/2.+a*ix-avx2;
	  for(int iy=0; iy<N; iy++)
	    {
	      pos = ix*N+iy;


	      if(param->getRunningCoupling())
		{
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2A = lat->cells[pos]->getg2mu2A();
		    }
		  else 
		    g2mu2A = 0;
		  
		  if(pos>0 && pos<(N-1)*N+N-1)
		    {
		      g2mu2B = lat->cells[pos]->getg2mu2B();
		    }
		  else
		    g2mu2B = 0;
		  
		  if(param->getRunWithQs()==2)
		    {
		      if(g2mu2A > g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==0)
		    {
		      if(g2mu2A < g2mu2B)
			Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		      else
			Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  else if(param->getRunWithQs()==1)
		    {
		      Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		    }
		  
		  if ( param->getRunWithLocalQs() == 1 )
		    {
		      // 3 flavors
		      alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		      //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		      
		      gfactor = g*g/(4.*PI*alphas);
		      //     cout << "Qs=" << Qs << endl;
		      //cout << alphas << endl;
		      // run with the local (in transverse plane) coupling
		    }
		  else
		    {
		      if ( param->getRunWithQs() == 0 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 1 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		      else if ( param->getRunWithQs() == 2 )
			alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		      
		      gfactor = g*g/(4.*PI*alphas);
		    }
		}
	      else
		gfactor = 1.;
	      
	      
 	      if (lat->cells[pos]->getEpsilon()*gfactor < cutoff)  //this is 1/fm^4, so Lambda_QCD^{-4}
		{
		  weight = 0.;
		}
	      else
		{
		  weight = lat->cells[pos]->getEpsilon()*gfactor;
		}
	      

	      y = -L/2.+a*iy-avy2;
	      avx += x*weight;
	      avy += y*weight;
	      avxSq += x*x*weight;
	      avySq += y*y*weight;
	      toteps += weight;
	    }
	}
 
      // cout << "toteps=" << toteps << endl;

      avx/=toteps;
      avy/=toteps;
      avxSq/=toteps;
      avySq/=toteps;
      avrSq/=toteps;
      Rbar = 1./sqrt(1./avxSq+1./avySq);
   
      //      cout << "new avx=" << avx << endl;
      //cout << "new avy=" << avy << endl;
      //cout << "sqrt(<r^2>)=" << sqrt(avrSq) << endl;
      //cout << "x_max=" << maxX << endl;
      //cout << "y_max=" << maxY << endl;
      //cout << "avxSq=" << avxSq << endl;
      //cout << "avySq=" << avySq << endl;
      //     cout << "Rbar=" << Rbar << " fm" << endl;

      param->setEccentricity2(eccentricity2);
      if (it==1)
	param->setPsi(Psi2);

      if(doAniso == 0)
	{
	  ofstream foutEcc(ecc_name.c_str(),ios::app); 
	  foutEcc << it*a*param->getdtau() << " " << eccentricity1 << " " << Psi1 << " " << eccentricity2 << " " << Psi2 << " " 
		  << eccentricity3 << " " << Psi3 << " " << eccentricity4 << " " << Psi4 
		  << " " << eccentricity5 << " " << Psi5 << " " << eccentricity6 << " " << Psi6 << " " << cutoff << " " << sqrt(avrSq) 
		  << " " << maxX << " " << maxY << " " << param->getb() << " " << param->getTpp() << " " << param->getArea() << " " << Rbar << endl;
	  foutEcc.close();      
	}

      if(doAniso==1)
	{
	  stringstream straniso_name;
	  straniso_name << "anisotropy" << param->getMPIRank() << ".dat";
	  string aniso_name;
	  aniso_name = straniso_name.str();
	  
	  ofstream foutAniso(aniso_name.c_str(),ios::app); 
	  
	  double TxxRot, TyyRot;
	  double ux,uy,PsiU;
	  double num=0., den=0.;
	  double unum=0., uden=0.;
	  double num2=0., den2=0.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  ux = lat->cells[pos]->getux();
		  uy = lat->cells[pos]->getuy();
		  unum += sqrt(ux*ux+uy*uy)*sin(2.*atan2(uy,ux));
		  uden += sqrt(ux*ux+uy*uy)*cos(2.*atan2(uy,ux));
		}
	    }

	  PsiU = atan2(unum,uden)/2.;

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  double Psi = PsiU;// param->getPsi();//-Pi/2.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << "Psi2=" << Psi2 << ", cos(Psi2)=" << cos(Psi2) << ", sin(Psi2)=" << sin(Psi2) << endl;
	  foutAniso << "PsiU=" << PsiU << ", cos(PsiU)=" << cos(PsiU) << ", sin(PsiU)=" << sin(PsiU) << endl;
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+Pi/8.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+Pi/4.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+3.*Pi/8.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+Pi/2.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 
	  
	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+5.*Pi/8.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+3.*Pi/4.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+7.*Pi/8.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+Pi;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 
	  
	  num=0.; den=0.;
	  num2=0.; den2=0.;
	  Psi = PsiU+9.*Pi/8.;

	  for(int ix=0; ix<N; ix++) 
	    {
	      for(int iy=0; iy<N; iy++)
		{
		  pos = (ix)*N+(iy);
		  
		  TxxRot = cos(Psi)*(cos(Psi)*lat->cells[pos]->getTxx()-sin(Psi)*lat->cells[pos]->getTxy())
		    -sin(Psi)*(cos(Psi)*lat->cells[pos]->getTxy()-sin(Psi)*lat->cells[pos]->getTyy());
		  TyyRot = sin(Psi)*(sin(Psi)*lat->cells[pos]->getTxx()+cos(Psi)*lat->cells[pos]->getTxy())
		    +cos(Psi)*(sin(Psi)*lat->cells[pos]->getTxy()+cos(Psi)*lat->cells[pos]->getTyy());
		  
		  num2 += lat->cells[pos]->getTxx()-lat->cells[pos]->getTyy();
		  den2 += lat->cells[pos]->getTxx()+lat->cells[pos]->getTyy();
		  
		  num += TxxRot-TyyRot;
		  den += TxxRot+TyyRot;
		}
	    }
	  
	  foutAniso << it*a*param->getdtau() << " " << num/den << " " << num2/den2 << " angle=" << Psi << endl; 

	  foutAniso.close();
	}
}


void Evolution::readNkt(Parameters *param)
{
  cout << "Reading n(k_T) from file ";
  string Npart, dummy;
  string kt, nkt, Tpp, b;
  double dkt;
  double dNdeta;

  // open file

  ifstream fin;
  stringstream strmult_name;
  strmult_name << "multiplicity" << param->getMPIRank() << ".dat";
  string mult_name;
  mult_name = strmult_name.str();
  fin.open(mult_name.c_str()); 
  cout << mult_name.c_str() << " ... " ;

  // open file

  ifstream fin2;
  stringstream strmult_name2;
  strmult_name2 << "NpartdNdy" << param->getMPIRank() << ".dat";
  string mult_name2;
  mult_name2 = strmult_name2.str();
  fin2.open(mult_name2.c_str()); 
  cout << mult_name2.c_str() << " ... " ;

  // read file

  if(fin)
    {
      for (int ikt=0; ikt<100; ikt++)
	{
	  if (!fin.eof())
	    {  
	      fin >> dummy;
	      fin >> kt;
	      fin >> nkt;
	      nIn[ikt]=atof(nkt.c_str());
	      fin >> dummy >> Tpp >> b >> Npart;
	      if(ikt==0)
		dkt = atof(kt.c_str());
	      if(ikt==1)
		dkt = dkt- atof(kt.c_str());
	    }
	  cout << nIn[ikt] << endl;
	}
      fin.close();
      cout << " done." << endl;
    }
  else
    {
      cout << "[Evolution.cpp:readNkt]: File " << mult_name.c_str() << " does not exist. Exiting." << endl;
      exit(1);
    }

  if(fin2)
    {
      if (!fin2.eof())
	{  
	  
	  fin2 >> Npart;
	  fin2 >> nkt;
	  fin2 >> Tpp;
	  fin2 >> b;
	  
	  dNdeta=atof(nkt.c_str());
	}
      fin2.close();
      cout << " done." << endl;
    }
  else
    {
      cout << "[Evolution.cpp:readNkt]: File " << mult_name2.c_str() << " does not exist. Exiting." << endl;
      exit(1);
    }


  double m,P;
  m=param->getJacobianm(); // in GeV
  P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
  double dNdeta2;
  dNdeta2 = 0.;

  for(int ik=0; ik<100; ik++)
    {
      if(param->getUsePseudoRapidity()==0)
	{
	  dNdeta2 += nIn[ik]*(ik+0.5)*dkt*dkt *2.*PI; //integrate, gives a ik*dkt*2pi*dkt
	}
      else
	{
	  dNdeta2 += nIn[ik]*(ik+0.5)*dkt*dkt *2.*PI
	    *cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt)*((ik+0.5)*dkt)))); 
	}
    }
  
  dNdeta *= cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/P/P));
  
  ofstream foutNN("NpartdNdy-mod.dat",ios::out); 
  foutNN << Npart << " " << dNdeta << " " << dNdeta2 << " " <<  atof(Tpp.c_str()) << " " <<  atof(b.c_str()) << endl;
  foutNN.close();


  exit(1);

}

  
int Evolution::multiplicity(Lattice *lat, Group *group, Parameters *param, int it)
{
  int N = param->getSize();
  int Nc = param->getNc();
  int npos, pos, posNew;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double Pi, kx, ky, kt2, omega2;
  double g = param->getg();
  Pi = param->getPi();
  int nn[2];
  nn[0] = N;
  nn[1] = N;
  double dtau = param->getdtau();
  double nkt;
  int bins = 100;
  double n[bins]; //k_T array
  double E[bins]; //k_T array
  double n2[bins]; //k_T array
  int counter[bins]; 
  double dkt = 2.83/static_cast<double>(bins);
  double dNdeta =0.;
  double dNdeta2 =0.;
  double dNdetaCut =0.;
  double dNdetaCut2 =0.;
  double dEdetaCut =0.;
  double dEdetaCut2 =0.;
  double dEdeta =0.;
  double dEdeta2 =0.;
  double Nkxky[N*N];
  
  stringstream strnkxky_name;
  strnkxky_name << "nkxky" << param->getMPIRank() << ".dat";
  string nkxky_name;
  nkxky_name = strnkxky_name.str();
 
  stringstream strNpartdNdy_name;
  strNpartdNdy_name << "NpartdNdy" << param->getMPIRank() << ".dat";
  string NpartdNdy_name;
  NpartdNdy_name = strNpartdNdy_name.str();

  stringstream strNpartdNdyH_name;
  strNpartdNdyH_name << "NpartdNdyHadrons" << param->getMPIRank() << ".dat";
  string NpartdNdyH_name;
  NpartdNdyH_name = strNpartdNdyH_name.str();

  stringstream strmult_name;
  strmult_name << "multiplicity" << param->getMPIRank() << ".dat";
  string mult_name;
  mult_name = strmult_name.str();

  stringstream strdNdy_name;
  strdNdy_name << "dNdy" << param->getMPIRank() << ".dat";
  string dNdy_name;
  dNdy_name = strdNdy_name.str();

  cout << "Measuring multiplicity ... " << endl;

  // fix transverse Coulomb gauge
  GaugeFix *gaugefix;
  gaugefix = new GaugeFix(nn);
  
  double maxtime;
  if ( param->getInverseQsForMaxTime() == 1 )
    {
      maxtime = 1./param->getAverageQs()*hbarc;
      cout << "maximal evolution time = " << maxtime << " fm" << endl; 
    }
  else
    {
      maxtime = param->getMaxtime(); // maxtime is in fm
    }

  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
  
  gaugefix->FFTChi(lat,group,param,4000);
  // gauge is fixed
  delete gaugefix;

  Matrix **E1;
  //  Matrix **E2;
  //  Matrix **pi;
  E1 = new Matrix*[N*N];
  //E2 = new Matrix*[N*N];
  //pi = new Matrix*[N*N];
  
  for(int i=0; i<N*N; i++)
    {
      E1[i] = new Matrix(Nc,0.);
      //  E2[i] = new Matrix(Nc,0.);
      // pi[i] = new Matrix(Nc,0.);
    }

  double g2mu2A, g2mu2B, gfactor, alphas, Qs;
  double c = param->getc();
  double muZero = param->getMuZero();

  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  
	  if(param->getRunningCoupling())
	    {
	  
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2A = lat->cells[pos]->getg2mu2A();
		}
	      else 
		g2mu2A = 0;
	      
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2B = lat->cells[pos]->getg2mu2B();
		}
	      else
		g2mu2B = 0;
	      
	      if(param->getRunWithQs()==2)
		{
		  if(g2mu2A > g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==0)
		{
		  if(g2mu2A < g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==1)
		{
		  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      
	      if ( param->getRunWithLocalQs() == 1)
		{	  
		  // 3 flavors
		  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		  //alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		  gfactor = g*g/(4.*PI*alphas);
		  // cout << "m Qs=" << Qs << endl;
		  // cout << alphas << endl;
		  // run with the local (in transverse plane) coupling
		}
	      else
		{
		  if ( param->getRunWithQs() == 0 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 1 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 2 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		  
		  gfactor = g*g/(4.*PI*alphas);
		}
	    }
	  else
	    gfactor = 1.;

	  if(param->getRunWithkt()==0)
	    {
	      *E1[pos] = lat->cells[pos]->getE1()*sqrt(gfactor); // replace one of the 1/g in the lattice E^i by the running one
	      //	      cout << "m " << sqrt(gfactor) << endl;

	    }
	  else
	    {
	      *E1[pos] = lat->cells[pos]->getE1(); 
	    }
	}
    }

  // do Fourier transforms
  fft->fftn(E1,E1,nn,2,1);
 
  for(int ik=0; ik<bins; ik++)
    {
      n[ik] = 0.;
      E[ik] = 0.;
      n2[ik] = 0.;
      counter[ik]=0;
    }
  
  int hbins = 2000;
  int kcounter=0;
  double Nhad = 0.;
  double Ehad = 0.;
  double Nh[hbins+1], Nhgsl[hbins+1], Eh[hbins+1], Ehgsl[hbins+1], NhL[hbins+1], NhLgsl[hbins+1], NhH[hbins+1], NhHgsl[hbins+1], Ng;
  for (int ih=0; ih<=hbins; ih++)
    {
      Nh[ih]=0.;
      Eh[ih]=0.;
      NhL[ih]=0.;
      NhH[ih]=0.;
    }

  // if(it==1)
  //   nkxky_name << "nkxkyInitial.dat";
  //  if (it==itmax)
  //  nkxky_name << "nkxkyFinal.dat";

  ofstream foutNkxky((nkxky_name).c_str(),ios::out); 
  
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
  	{
  	  pos = i*N+j;
  	  npos = (N-i)*N+(N-j);
	  
  	  kx = 2.*param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N));
  	  ky = 2.*param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N));

  	  //	  if (sqrt(kx*kx+ky*ky)>2.)
  	  //continue;

  	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));//
  	  omega2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice dispersion relation (this is omega squared)

  	  // i=0 or j=0 have no negative k_T value available
 
  	  if(i!=0 && j!=0)
  	    {
  	      //cout << ((*E1[pos])*(*E1[npos])).trace()  << endl;
  	      //cout << ((*E2[pos])*(*E2[npos])).trace()  << endl;
  	      //cout << ((*pi[pos])*(*pi[npos])).trace()  << endl << endl;

  	      if(omega2!=0)
  		{
  		  nkt = 2./sqrt(omega2)/static_cast<double>(N*N) * ( g*g/((it-0.5)*dtau)*( (((*E1[pos])*(*E1[npos])).trace()).real() 
											   )
  								     );
  		  if(param->getRunWithkt()==1)
  		    {
  		      nkt *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
  		    }
  		}
	      
	      
  	      dNdeta += nkt;
  	      dEdeta += nkt*sqrt(omega2)*0.1973269718/a;
	          
  	      for(int ik=0; ik<bins; ik++)
  		{
  		  if (abs(sqrt(kt2))>ik*dkt && abs(sqrt(kt2))<=(ik+1)*dkt)
  		    {
  		      n[ik]+=nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      E[ik]+=sqrt(omega2)*0.1973269718/a *nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      n2[ik]+=nkt/dkt/2/Pi/sqrt(kt2);
  		      // dividing by bin size; bin is dkt times Jacobian k(=ik*dkt) times 2Pi in phi 
  		      // times the correct number of counts for an infinite lattice: area in bin divided by total area
  		      counter[ik]+=1; // number of entries in n[ik]
  		    }
  		}
  	    }
	  Nkxky[pos] = nkt*N*N/Pi/Pi/2./2.;
  	}
    }
  //foutNkxky.close();

  
  

  /// -------- 2 ---------

  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  
	  if(param->getRunningCoupling())
	    {
	  
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2A = lat->cells[pos]->getg2mu2A();
		}
	      else 
		g2mu2A = 0;
	      
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2B = lat->cells[pos]->getg2mu2B();
		}
	      else
		g2mu2B = 0;
	      
	      if(param->getRunWithQs()==2)
		{
		  if(g2mu2A > g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==0)
		{
		  if(g2mu2A < g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==1)
		{
		  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      
	      if ( param->getRunWithLocalQs() == 1)
		{	  
		  // 3 flavors
		  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		  //alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		  gfactor = g*g/(4.*PI*alphas);
		  // cout << "m Qs=" << Qs << endl;
		  // cout << alphas << endl;
		  // run with the local (in transverse plane) coupling
		}
	      else
		{
		  if ( param->getRunWithQs() == 0 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 1 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 2 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		  
		  gfactor = g*g/(4.*PI*alphas);
		}
	    }
	  else
	    gfactor = 1.;

	  if(param->getRunWithkt()==0)
	    {
	      *E1[pos] = lat->cells[pos]->getE2()*sqrt(gfactor); // "
	    }
	  else
	    {
	      *E1[pos] = lat->cells[pos]->getE2(); 
	    }
	}
    }

  fft->fftn(E1,E1,nn,2,1);
  
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
  	{
  	  pos = i*N+j;
  	  npos = (N-i)*N+(N-j);
	  
  	  kx = 2.*param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N));
  	  ky = 2.*param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N));

  	  //	  if (sqrt(kx*kx+ky*ky)>2.)
  	  //continue;

  	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));//
  	  omega2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice dispersion relation (this is omega squared)

  	  // i=0 or j=0 have no negative k_T value available
 
  	  if(i!=0 && j!=0)
  	    {
  	      //cout << ((*E1[pos])*(*E1[npos])).trace()  << endl;
  	      //cout << ((*E2[pos])*(*E2[npos])).trace()  << endl;
  	      //cout << ((*pi[pos])*(*pi[npos])).trace()  << endl << endl;

  	      if(omega2!=0)
  		{
  		  nkt = 2./sqrt(omega2)/static_cast<double>(N*N) * ( g*g/((it-0.5)*dtau)*( ((((*E1[pos])*(*E1[npos])).trace()).real())));
		  if(param->getRunWithkt()==1)
  		    {
  		      nkt *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
  		    }
  		}
	      
  	      dNdeta += nkt;
  	      dEdeta += nkt*sqrt(omega2)*0.1973269718/a;
	          
  	      for(int ik=0; ik<bins; ik++)
  		{
  		  if (abs(sqrt(kt2))>ik*dkt && abs(sqrt(kt2))<=(ik+1)*dkt)
  		    {
  		      n[ik]+=nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      E[ik]+=sqrt(omega2)*0.1973269718/a *nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      n2[ik]+=nkt/dkt/2/Pi/sqrt(kt2);
  		      // dividing by bin size; bin is dkt times Jacobian k(=ik*dkt) times 2Pi in phi 
  		      // times the correct number of counts for an infinite lattice: area in bin divided by total area
		      // counter[ik]+=1; // number of entries in n[ik]
  		    }
  		}
  	    }
	  Nkxky[pos] += nkt*N*N/Pi/Pi/2./2.;
	  //	  if (param->getWriteOutputs()==2)
	  //   foutNkxky << i << " " << j << " " << nkt << "\n";
  	}
      //      if (param->getWriteOutputs()==2)
      // foutNkxky << endl;
    }
  //foutNkxky.close();
  

  /// ------3 --------

  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  
	  if(param->getRunningCoupling())
	    {
	  
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2A = lat->cells[pos]->getg2mu2A();
		}
	      else 
		g2mu2A = 0;
	      
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2B = lat->cells[pos]->getg2mu2B();
		}
	      else
		g2mu2B = 0;
	      
	      if(param->getRunWithQs()==2)
		{
		  if(g2mu2A > g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==0)
		{
		  if(g2mu2A < g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==1)
		{
		  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      
	      if ( param->getRunWithLocalQs() == 1)
		{	  
		  // 3 flavors
		  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		  //alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		  gfactor = g*g/(4.*PI*alphas);
		  // cout << "m Qs=" << Qs << endl;
		  // cout << alphas << endl;
		  // run with the local (in transverse plane) coupling
		}
	      else
		{
		  if ( param->getRunWithQs() == 0 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 1 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 2 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		  
		  gfactor = g*g/(4.*PI*alphas);
		}
	    }
	  else
	    gfactor = 1.;

	  if(param->getRunWithkt()==0)
	    {
	      *E1[pos] = lat->cells[pos]->getpi()*sqrt(gfactor); // replace the only 1/g by the running one (physical pi goes like 1/g, like physical E^i)
	      //	      cout << "m " << sqrt(gfactor) << endl;

	    }
	  else
	    {
	      *E1[pos] = lat->cells[pos]->getpi(); 
	    }
	}
    }

  // do Fourier transforms
  fft->fftn(E1,E1,nn,2,1);

  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
  	{
  	  pos = i*N+j;
  	  npos = (N-i)*N+(N-j);
	  
  	  kx = 2.*param->getPi()*(-0.5+static_cast<double>(i)/static_cast<double>(N));
  	  ky = 2.*param->getPi()*(-0.5+static_cast<double>(j)/static_cast<double>(N));

  	  //	  if (sqrt(kx*kx+ky*ky)>2.)
  	  //continue;

  	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));//
  	  omega2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice dispersion relation (this is omega squared)

  	  // i=0 or j=0 have no negative k_T value available
 
  	  if(i!=0 && j!=0)
  	    {
  	      //cout << ((*E1[pos])*(*E1[npos])).trace()  << endl;
  	      //cout << ((*E2[pos])*(*E2[npos])).trace()  << endl;
  	      //cout << ((*pi[pos])*(*pi[npos])).trace()  << endl << endl;

  	      if(omega2!=0)
  		{
  		  nkt = 2./sqrt(omega2)/static_cast<double>(N*N) * ( ((it-0.5)*dtau) * ( (((*E1[pos])*(*E1[npos])).trace()).real() ) );
  		  if(param->getRunWithkt()==1)
  		    {
  		      nkt *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
  		    }
  		}
	      
	      
  	      dNdeta += nkt;
  	      dEdeta += nkt*sqrt(omega2)*0.1973269718/a;
	          
  	      for(int ik=0; ik<bins; ik++)
  		{
  		  if (abs(sqrt(kt2))>ik*dkt && abs(sqrt(kt2))<=(ik+1)*dkt)
  		    {
  		      n[ik]+=nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      E[ik]+=sqrt(omega2)*0.1973269718/a *nkt/dkt/2/Pi/sqrt(kt2) * 2*Pi*sqrt(kt2)*dkt *N*N/Pi/Pi/2./2.;
  		      n2[ik]+=nkt/dkt/2/Pi/sqrt(kt2);
  		      // dividing by bin size; bin is dkt times Jacobian k(=ik*dkt) times 2Pi in phi 
  		      // times the correct number of counts for an infinite lattice: area in bin divided by total area
  		      //counter[ik]+=1; // number of entries in n[ik]
  		    }
  		}
  	    }
	  Nkxky[pos] += nkt *N*N/Pi/Pi/2./2.;
	  if (param->getWriteOutputs()==2)
	    foutNkxky << 2.*sin(kx/2.)/a*0.1973269718 << " " << 2.*sin(ky/2.)/a*0.1973269718 << " " << Nkxky[pos]*a/0.1973269718*a/0.1973269718 << "\n";
  	}
      if (param->getWriteOutputs()==2)
	foutNkxky << endl;
    }
  foutNkxky.close();
  




  double m,P;
  m=param->getJacobianm(); // in GeV
  P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV

  ofstream foutMult(mult_name.c_str(),ios::app); 
  //  ofstream foutdNdy(dNdy_name.c_str(),ios::app); 
  for(int ik=0; ik<bins; ik++)
    {
      if(counter[ik]>0)
	{
	  n[ik] = n[ik]/static_cast<double>(counter[ik]);
	  E[ik] = E[ik]/static_cast<double>(counter[ik]);
	  if(param->getUsePseudoRapidity()==0)
	    {
	      dNdeta2 += n[ik]*(ik+0.5)*dkt*dkt *2.*PI; //integrate, gives a ik*dkt*2pi*dkt
	      dEdeta2 += E[ik]*(ik+0.5)*dkt*dkt *2.*PI; //integrate, gives a ik*dkt*2pi*dkt
	      if(ik*dkt/a*0.1973269718>3.) // 
		{
		  dNdetaCut += n[ik]*(ik+0.5)*dkt*dkt *2.*PI; 
		  dEdetaCut += E[ik]*(ik+0.5)*dkt*dkt *2.*PI; 
		}
	      if(ik*dkt/a*0.1973269718>6.) // large cut
		{
		  dNdetaCut2 += n[ik]*(ik+0.5)*dkt*dkt *2.*PI; 
		  dEdetaCut2 += E[ik]*(ik+0.5)*dkt*dkt *2.*PI; 
		}
	    }
	  else
	    {
	      dNdeta2 += n[ik]*(ik+0.5)*dkt*dkt *2.*PI
		*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718))));  //integrate, gives a ik*dkt*2pi*dkt
	      dEdeta2 += E[ik]*(ik+0.5)*dkt*dkt *2.*PI*
		cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718))));  //integrate, gives a ik*dkt*2pi*dkt

	      if(ik*dkt/a*0.1973269718>3.) //
		{
		  dNdetaCut += n[ik]*(ik+0.5)*dkt*dkt *2.*PI 
		    *cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718))));     
		    dEdetaCut += E[ik]*(ik+0.5)*dkt*dkt *2.*PI
		    *cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718)))); 
		}
	      if(ik*dkt/a*0.1973269718>6.) // large cut
		{
		  dNdetaCut2 += n[ik]*(ik+0.5)*dkt*dkt *2.*PI
		    *cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718)))); 
		  dEdetaCut2 += E[ik]*(ik+0.5)*dkt*dkt *2.*PI
		    *cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(((ik+0.5)*dkt/a*0.1973269718)*((ik+0.5)*dkt/a*0.1973269718)))); 
		}
	  
	    }
	  //integrate, gives a ik*dkt*2pi*dkt, in |eta|<2.4, 0.4 GeV p_T cut, charged N_track (offline, factor 0.83)
	}
      
      // output dN/d^2k
      if(it == itmax)
	{
	  foutMult << it*dtau*a << " " << ik*dkt/a*0.1973269718 << " " 
		   << n[ik]*a/0.1973269718*a/0.1973269718 << " " 
		   << n2[ik]*a/0.1973269718*a/0.1973269718 << " " << param->getTpp() << " " << param->getb() << " " << param->getNpart() << endl;
	}
    }
  
  //*ik*dkt*ik*dkt*ik*dkt*ik*dkt/a*0.1973269718/a*0.1973269718/a*0.1973269718/a*0.1973269718 // this is a factor of k_T^4

  //  foutdNdy << it*dtau*a << " " << dNdeta << " " << dNdeta2 << " " << dNdetaCut << " " << dEdeta << " " << dEdeta2 << endl;
  // in this version the gfactor is included above (and depends on the position)
  //foutdNdy.close();
  foutMult.close();

  
  double Ech=0.;
  double Ech2=0.;
  double Nch=0.;
  int icheck=0;
  double dNdetaHadrons, dNdetaHadronsCut, dNdetaHadronsCut2;
  double dEdetaHadrons, dEdetaHadronsCut, dEdetaHadronsCut2;

  // compute hadrons using fragmentation function 
  if(it == itmax)
    {
      //      cout << " before Hadronizing ... " << endl;
      //sleep(10);
      cout << " Hadronizing ... " << endl;
      double z,frac;
      double mypt, kt;
      int ik;
      int steps=6000;
      double dz = 0.95/static_cast<double>(steps);
      double zValues[steps+1];
      double zintegrand[steps+1];
      //double Ezintegrand[steps+1];
      //double Lzintegrand[steps+1];
      //double Hzintegrand[steps+1];
      gsl_interp_accel *zacc = gsl_interp_accel_alloc ();
      gsl_spline *zspline = gsl_spline_alloc (gsl_interp_cspline, steps+1);

      for(int ih=0; ih<=hbins; ih++)
	{

	  mypt = ih * (20./static_cast<double>(hbins));  
	  
	  for (int iz=0; iz<=steps; iz++)
	    {
	      z = 0.05 + iz * dz;
	      zValues[iz] = z;
	      
	      kt = mypt/z;

	      ik = static_cast<int>(floor(kt*a/0.1973269718/dkt-0.5+0.00000001));
	      
	      frac = (kt - (ik+0.5)*dkt/a*0.1973269718)/(dkt/a*0.1973269718);
     
	      if(ik+1<bins && ik >=0)
		Ng = ((1.-frac)*n[ik]+frac*n[ik+1])*a/0.1973269718*a/0.1973269718; // to make dN/d^2k_T fo k_T in GeV 
	      else 
		Ng =0.;

	      
	      if(param->getUsePseudoRapidity()==0)
		{
		  zintegrand[iz] = 1./(z*z) * Ng * frag->kkp(7,1,z,kt);
		  // Ezintegrand[iz] = mypt * 1./(z*z) * Ng * frag->kkp(7,1,z,kt);
		  // Lzintegrand[iz] = 1./(z*z) * Ng * frag->kkp(7,1,z,kt/2.);
		  // Hzintegrand[iz] = 1./(z*z) * Ng * frag->kkp(7,1,z,kt*2.);
		}
	      else
		{
		  zintegrand[iz] =  1./(z*z) * Ng * 
			2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
			      +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
			      +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))));
	
		  // Ezintegrand[iz] =  mypt * 1./(z*z) * Ng * 
		  // 	2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
		  // 	      +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
		  // 	      +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))));	
		  
		  // Lzintegrand[iz] =  1./(z*z) * Ng * 
		  // 	2. * (frag->kkp(1,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
		  // 	      +frag->kkp(2,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
		  // 	      +frag->kkp(4,1,z,kt/2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))));

		  // Hzintegrand[iz] =  1./(z*z) * Ng * 
		  // 	2. * (frag->kkp(1,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
		  // 	      +frag->kkp(2,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
		  // 	      +frag->kkp(4,1,z,kt*2.)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))));

		  
			}
	    }

	  zValues[steps]=1.; //set exactly 1
	  
	  
	  gsl_spline_init (zspline, zValues, zintegrand, steps+1);
	  Nhgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
	  

	  //	  gsl_spline_init (zspline, zValues, Lzintegrand, steps+1);
	  //NhLgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);
	  
	  //gsl_spline_init (zspline, zValues, Hzintegrand, steps+1);
	  //NhHgsl[ih] = gsl_spline_eval_integ(zspline, 0.05, 1., zacc);

	}

      gsl_spline_free(zspline);
      gsl_interp_accel_free(zacc);
          
      stringstream strmultHad_name;
      strmultHad_name << "multiplicityHadrons" << param->getMPIRank() << ".dat";
      string multHad_name;
      multHad_name = strmultHad_name.str();
      
      ofstream foutdNdpt(multHad_name.c_str(),ios::out); 
      for (int ih=0; ih<=hbins; ih++)
	{
	  if (ih%10==0)
	foutdNdpt << ih * 20. /static_cast<double>(hbins)  << " " << Nhgsl[ih] << " " << 0. << " " << 0. 
		  << " " << param->getTpp() << " " << param->getb() << endl; // leaving out the L and H ones for now
	}
      foutdNdpt.close();
      
      cout << " done." << endl;
    
      //integrate over pT using gsl      
      double pt[hbins+1];
      double integrand[hbins+1];
      double Eintegrand[hbins+1];
      for(int ih=0; ih<=hbins; ih++)
	{
	  pt[ih] = ih * 20. /static_cast<double>(hbins);
	  integrand[ih] = Nhgsl[ih]*pt[ih];
	  Eintegrand[ih] = Nhgsl[ih]*pt[ih]*pt[ih];
	}

      gsl_interp_accel *ptacc = gsl_interp_accel_alloc ();
      gsl_spline *ptspline = gsl_spline_alloc (gsl_interp_cspline, hbins+1);
      gsl_spline_init (ptspline, pt, integrand, hbins+1);
      dNdetaHadrons = 2*Pi*gsl_spline_eval_integ(ptspline, 0.25, 19., ptacc);
      dNdetaHadronsCut = 2*Pi*gsl_spline_eval_integ(ptspline, 3., 19., ptacc);
      dNdetaHadronsCut2 = 2*Pi*gsl_spline_eval_integ(ptspline, 6., 19., ptacc);

      gsl_spline_init (ptspline, pt, Eintegrand, hbins+1);
      dEdetaHadrons = 2*Pi*gsl_spline_eval_integ(ptspline, 0.25, 19., ptacc);
      dEdetaHadronsCut = 2*Pi*gsl_spline_eval_integ(ptspline, 3., 19., ptacc);
      dEdetaHadronsCut2 = 2*Pi*gsl_spline_eval_integ(ptspline, 6., 19., ptacc);

      gsl_spline_free(ptspline);
      gsl_interp_accel_free(ptacc);
      //      cout << "gsl integral: dNdeta= " << dNdetaHadrons << " dEdeta= " << dEdetaHadrons << endl;
      
      // icheck=0.;
      // // compute N_ch and mean p_T
      // for (int ih=0; ih<=hbins; ih++)
      // 	{
      // 	  if ( ih * 20. /static_cast<double>(hbins) > 0.25 &&  ih * 20. /static_cast<double>(hbins) < 20.) // p_T cut in GeV
      // 	    {
      // 	      if ( icheck==0 || ih==hbins )
      // 		{
      // 		  Nch += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) * Nh[ih] * (20. /static_cast<double>(hbins)) * 0.5;
      // 		  Ech += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) *
      // 		    Nh[ih] * ( ih * 20. /static_cast<double>(hbins) ) * (20. /static_cast<double>(hbins)) * 0.5;
      // 		  Ech2 += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) * Eh[ih]  * (20. /static_cast<double>(hbins)) * 0.5;
      // 		  icheck=1;
      // 		}
      // 	      else
      // 		{
      // 		  Nch += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) * Nh[ih] * (20. /static_cast<double>(hbins));
      // 		  Ech += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) * 
      // 		    Nh[ih] * ( ih * 20. /static_cast<double>(hbins) ) * (20. /static_cast<double>(hbins)) ;
      // 		  Ech2 += 2*Pi* ( ih * 20. /static_cast<double>(hbins) ) * Eh[ih]  * (20. /static_cast<double>(hbins));
      // 		}
      // 	    }
      // 	}  
      
      // cout << "looser intergral: dNdeta= " << Nch << " dEdeta= " << Ech << endl;
    }
  
  if(param->getUsePseudoRapidity()==0 && param->getMPIRank()==0)
    {
      cout << "dN/dy 1 = " << dNdeta << ", dE/dy 1 = " << dEdeta << endl; 
	  cout << "dN/dy 2 = " << dNdeta2 << ", dE/dy 2 = " << dEdeta2 << endl; 
	  cout << "gluon <p_T> = " << dEdeta/dNdeta << endl;
    }
  else
    {
      m=param->getJacobianm(); // in GeV
      P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
      dNdeta*=cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P)));
      dEdeta*=cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+m*m/(P*P)));
      
      if(param->getMPIRank()==0)
	{
	  cout << "dN/deta 1 = " << dNdeta << ", dE/deta 1 = " << dEdeta << endl; 
	  cout << "dN/deta 2 = " << dNdeta2 << ", dE/deta 2 = " << dEdeta2 << endl; 
	  cout << "dN/deta_cut 1 = " << dNdetaCut << endl; 
	  cout << "dN/deta_cut 2 = " << dNdetaCut2 << endl; 
	  cout << "gluon <p_T> = " << dEdeta/dNdeta << endl;
	}
    }

  if (dNdeta == 0.)
    {
      cout << "No collision happened on rank " << param->getMPIRank() << ". Restarting with new random number..." << endl;
      //MPI::Finalize();
      //exit(1);
      for(int i=0; i<N*N; i++)
	{
	  delete E1[i];
	  //	  delete E2[i];
	  // delete pi[i];
	}
      
      delete[] E1;
      //      delete[] E2;
      // delete[] pi;
      
      //      delete gaugefix;
      return 0;
    }
  

  if (it==itmax)
    {
      cout << "hadron <p_T> = " << dEdetaHadrons/dNdetaHadrons << endl;
      cout << "Hadrons: dN/dy(p_T>250 MeV)=" << dNdetaHadrons << ", dE/dy(p_T>250 MeV)=" << dEdetaHadrons << endl;
      
      stringstream strmeanpt_name;
      strmeanpt_name << "meanpt" << param->getMPIRank() << ".dat";
      string meanpt_name;
      meanpt_name = strmeanpt_name.str();

      ofstream foutNch(meanpt_name.c_str(),ios::out); 
      foutNch << dNdeta << " " << dEdeta/dNdeta << " " << dNdetaHadrons << " " << dEdetaHadrons/dNdetaHadrons << endl;
      foutNch.close();
      
      ofstream foutNN(NpartdNdy_name.c_str(),ios::app); 
      foutNN << param->getNpart() << " " << dNdeta << " " << param->getTpp() << " " << param->getb() << " " << dEdeta << " " << param->getRandomSeed() 
	     << " " <<  "N/A" << " " << "N/A" << " " << "N/A" << " " <<  dNdetaCut << " " << dEdetaCut 
	     << " " << dNdetaCut2 << " " << dEdetaCut2 << " " << g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)))) << endl;
      foutNN.close();

      ofstream foutNNH(NpartdNdyH_name.c_str(),ios::app); 
      foutNNH << param->getNpart() << " " << dNdetaHadrons << " " << param->getTpp() << " " << param->getb() << " " << dEdetaHadrons << " " << param->getRandomSeed()   << " " <<  "N/A" << " " << "N/A" << " " << "N/A" << " " <<  dNdetaHadronsCut << " " << dEdetaHadronsCut 
	      << " " << dNdetaHadronsCut2 << " " << dEdetaHadronsCut2 << " " << g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)))) << endl;
      foutNNH.close();

    }
  
  for(int i=0; i<N*N; i++)
    {
      delete E1[i];
      //     delete E2[i];
      // delete pi[i];
    }
  
  delete[] E1;
  //  delete[] E2;
  //delete[] pi;

  //  delete gaugefix;
  cout << " done." << endl;
  param->setSuccess(1);
  return 1;

}



int Evolution::correlations(Lattice *lat, Group *group, Parameters *param, int it)
{
  int N = param->getSize();
  int Nc = param->getNc();
  int npos, pos, posNew;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double Pi, kx, ky, kt2, omega2;
  double g = param->getg();
  Pi = param->getPi();
  int nn[2];
  nn[0] = N;
  nn[1] = N;
  double dtau = param->getdtau();
  double nkt, nkt1, nkt2, nkt3, nkt4, nkt5, nkt6, nktNoMixedTerms;
  int bins = 40;
  int phiBins = 16;
  double n[bins][phiBins]; // |k_T|, phi array
  double n2[bins][phiBins]; // |k_T|, phi array
  double nkxky[N][N]; // kx, ky array
  double nk[bins]; //|k_T| array
  double nNoMixedTerms[bins][phiBins]; //|k_T|, phi array
  double nkNoMixedTerms[bins]; //|k_T| array
  int counter[bins][phiBins]; 
  int counterk[bins]; 
  double dkt = 2.83/static_cast<double>(bins);
  double dNdeta =0.;
  double dNdetaNoMixedTerms =0.;
  double dNdeta1 =0.;
  double dNdeta2 =0.;
  double dNdeta3 =0.;
  double dNdeta4 =0.;
  double dNdeta5 =0.;
  double dNdeta6 =0.;
  double anglePhi;
  double k;  
  double deltaPhi = 2.*Pi/static_cast<double>(phiBins);

  stringstream strCorr_name;
  strCorr_name << "Corr" << param->getMPIRank() << ".dat";
  string Corr_name;
  Corr_name = strCorr_name.str();

  stringstream strPhiMult_name;
  strPhiMult_name << "MultPhi" << param->getMPIRank() << ".dat";
  string PhiMult_name;
  PhiMult_name = strPhiMult_name.str();

  stringstream strPhiMultHad_name;
  strPhiMultHad_name << "MultPhiHadrons" << param->getMPIRank() << ".dat";
  string PhiMultHad_name;
  PhiMultHad_name = strPhiMultHad_name.str();

  cout << "Measuring multiplicity version 2... " << endl;

  // fix transverse Coulomb gauge
  GaugeFix *gaugefix;
  gaugefix = new GaugeFix(nn);
  
  double maxtime;
  if ( param->getInverseQsForMaxTime() == 1 )
    {
      maxtime = 1./param->getAverageQs()*hbarc;
      cout << "maximal evolution time = " << maxtime << " fm" << endl; 
    }
  else
    {
      maxtime = param->getMaxtime(); // maxtime is in fm
    }

  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
  gaugefix->FFTChi(lat,group,param,4000);
 
  // gauge is fixed
  Matrix U1(Nc,1.);
  Matrix U2(Nc,1.);
  Matrix U1dag(Nc,1.);
  Matrix U2dag(Nc,1.);


 
  Matrix **A1;
  A1 = new Matrix*[N*N];
  Matrix **A2;
  A2 = new Matrix*[N*N];
  Matrix **phi;
  phi = new Matrix*[N*N];
  
  Matrix **E1;
  Matrix **E2;
  Matrix **pi;
  E1 = new Matrix*[N*N];
  E2 = new Matrix*[N*N];
  pi = new Matrix*[N*N];
  
  for(int i=0; i<N*N; i++)
    {
      A1[i] = new Matrix(Nc,0.);
      A2[i] = new Matrix(Nc,0.);
      E1[i] = new Matrix(Nc,0.);
      E2[i] = new Matrix(Nc,0.);
      pi[i] = new Matrix(Nc,0.);
      phi[i] = new Matrix(Nc,0.);
    }
  
  
 //  for (int i=0; i<N; i++)
//     {
//       for (int j=0; j<N; j++)
// 	{
// 	  pos = i*N+j;
// 	  U1 = lat->cells[pos]->getUx();
// 	  U2 = lat->cells[pos]->getUy();

// 	  U1dag = U1;
// 	  U2dag = U2;
// 	  U1dag.conjg();
// 	  U2dag.conjg();
	
// 	  *A1[pos] = complex<double>(0.,-0.5)*(U1-U1dag);
// 	  *A2[pos] = complex<double>(0.,-0.5)*(U2-U2dag);

// 	  if (i==N/2 && j==N/2)
// 	    cout << "approximation=" << endl << *A1[pos] << endl << endl;
	  
// 	}
//     }

 //version that determines the exact log of U1 and U2:
 for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  U1 = lat->cells[pos]->getUx();
	  U2 = lat->cells[pos]->getUy();

	  //	  if( i==N/2 && j==N/2)
	  //       cout << "U1=" << U1 << endl;
	  
	  U1.logm();
	  U2.logm();

	//   if( i==N/2 && j==N/2)
// 	    {
// 	      cout << "log(U1) = " << U1 << endl;    
// 	      cout << "U1(retrieved)=" << U1.expm() << endl;
// 	      U1.logm();
// 	    }

	  *A1[pos] = complex<double>(0.,-1.)*U1;
	  *A2[pos] = complex<double>(0.,-1.)*U2;

// 	  if (i==N/2 && j==N/2)
// 	    cout << "exact=" << endl << *A1[pos] << endl << endl;

	}
    }



  double g2mu2A, g2mu2B, gfactor, alphas, Qs;
  double c = param->getc();
  double muZero = param->getMuZero();

  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  
	  if(param->getRunningCoupling())
	    {
	  
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2A = lat->cells[pos]->getg2mu2A();
		}
	      else 
		g2mu2A = 0;
	      
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2B = lat->cells[pos]->getg2mu2B();
		}
	      else
		g2mu2B = 0;
	      
	      if(param->getRunWithQs()==2)
		{
		  if(g2mu2A > g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==0)
		{
		  if(g2mu2A < g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==1)
		{
		  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      
	      if ( param->getRunWithLocalQs() == 1)
		{
		  // 3 flavors
		  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		  //alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
		  gfactor = g*g/(4.*PI*alphas);
		  //     cout << "Qs=" << Qs << endl;
		  //cout << alphas << endl;
		  // run with the local (in transverse plane) coupling
		}
	      else
		{
		  if ( param->getRunWithQs() == 0 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 1 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 2 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		  
		  gfactor = g*g/(4.*PI*alphas);
		}
	    }
	  else
	    gfactor = 1.;

	  if(param->getRunWithkt()==0)
	    {
	      *E1[pos] = lat->cells[pos]->getE1()*sqrt(gfactor); // replace one of the 1/g in the lattice E^i by the running one
	      *E2[pos] = lat->cells[pos]->getE2()*sqrt(gfactor); // "
	      *pi[pos] = lat->cells[pos]->getpi()*sqrt(gfactor); // replace the only 1/g by the running one (physical pi goes like 1/g, like physical E^i)
	      *A1[pos] = *A1[pos]*sqrt(gfactor); //"
	      *A2[pos] = *A2[pos]*sqrt(gfactor); // "
	      *phi[pos] = lat->cells[pos]->getphi()*sqrt(gfactor); // replace the only 1/g by the running one (physical pi goes like 1/g, like physical E^i)
	    }
	  else
	    {
	      *E1[pos] = lat->cells[pos]->getE1(); 
	      *E2[pos] = lat->cells[pos]->getE2(); 
	      *pi[pos] = lat->cells[pos]->getpi(); 
	      *phi[pos] = lat->cells[pos]->getphi(); 
	    }
	}
    }



  // do Fourier transforms

  fft->fftn(A1,A1,nn,2,1);
  fft->fftn(A2,A2,nn,2,1);
  fft->fftn(phi,phi,nn,2,1);

  fft->fftn(E1,E1,nn,2,1);
  fft->fftn(E2,E2,nn,2,1);
  fft->fftn(pi,pi,nn,2,1);


//   double x,y;
//     if(it==itmax)
// 	{	  
// 	  stringstream strepsx_name;
// 	  strepsx_name << "eps-x" << param->getMPIRank() << ".dat";
// 	  string epsx_name;
// 	  epsx_name = strepsx_name.str();

// 	  ofstream foutEpsX(epsx_name.c_str(),ios::app); 
// 	  for(int ix=1; ix<N; ix++) // loop over all positions
// 	    {
// 	      for(int iy=1; iy<N; iy++)
// 		{
// 		  pos = ix*N+iy;
// 		  x = -L/2.+a*ix;
// 		  y = -L/2.+a*iy;
		  
// 		  if(param->getRunningCoupling())
// 		    {
// 		      if(pos>0 && pos<(N-1)*N+N-1)
// 			{
// 			  g2mu2A = lat->cells[pos]->getg2mu2A();
// 			}
// 		      else 
// 			g2mu2A = 0;
		      
// 		      if(pos>0 && pos<(N-1)*N+N-1)
// 			{
// 			  g2mu2B = lat->cells[pos]->getg2mu2B();
// 			}
// 		      else
// 			g2mu2B = 0;
		      
// 		      if(param->getRunWithQs()==2)
// 			{
// 			  if(g2mu2A > g2mu2B)
// 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			  else
// 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			}
// 		      else if(param->getRunWithQs()==0)
// 			{
// 			  if(g2mu2A < g2mu2B)
// 			    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			  else
// 			    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			}
// 		      else if(param->getRunWithQs()==1)
// 			{
// 			  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
// 			}
		      

// 		      if ( param->getRunWithLocalQs() == 1) 
// 			{
// 			  // 3 flavors
// 			  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
// 			  //	      alphas = 12.*PI/((27.)*2.*log(param->getRunWithThisFactorTimesQs()*Qs/0.2)); // 3 flavors
			  
// 			  gfactor = g*g/(4.*PI*alphas);
// 			  //     cout << "Qs=" << Qs << endl;
// 			  //cout << alphas << endl;
// 			  // run with the local (in transverse plane) coupling
// 			}
// 		      else
// 			{
// 			  if ( param->getRunWithQs() == 0 )
// 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
// 			  else if ( param->getRunWithQs() == 1 )
// 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
// 			  else if ( param->getRunWithQs() == 2 )
// 			    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
			  
// 			  gfactor = g*g/(4.*PI*alphas);
// 			}
// 		    }
//                  else
// 		   gfactor = 1.;
  
// 		  foutEpsX << x << " " << y << " " << 1./static_cast<double>(N*N) * (complex<double>(0.,1.)*((*pi[pos])*(*phi[npos])-(*phi[pos])*(*pi[npos])).trace()).real()
// 			   << " " <<  1./static_cast<double>(N*N) * (complex<double>(0.,1.)*((*E1[pos])*(*A1[npos])-(*A1[pos])*(*E1[npos])+
// 										(*E2[pos])*(*A2[npos])-(*A2[pos])*(*E2[npos])).trace()).real()
// 			   << endl;   
//                  // abs just to get rid of negative 10^(-17) numbers at edge
// 		}
// 	      foutEpsX << endl;
// 	    }
// 	  foutEpsX << endl;
// 	  foutEpsX.close();
//  	}



  for(int ik=0; ik<bins; ik++)
    {
      nk[ik] = 0.;
      nkNoMixedTerms[ik] = 0.;
      counterk[ik] = 0;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  n[ik][iphi] = 0.;
	  n2[ik][iphi] = 0.;
	  nNoMixedTerms[ik][iphi] = 0.;
	  counter[ik][iphi]=0;
	}
    }

  // leave out the first cell to make it symmetric 
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  npos = (N-i)*N+(N-j);
	  
	  kx = 2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(N));
	  ky = 2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(N));
	  
	  // don't use momenta that are badly described on the lattice
	  //	  if (sqrt(kx*kx+ky*ky)>2.)
	  // continue;
	  
	  //kt2 = kx*kx+ky*ky;//
	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));//
	  //omega2 = kx*kx+ky*ky;
	  omega2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice dispersion relation (this is omega squared)
	  
	  // i=0 or j=0 have no negative k_T value available

	  if(i!=0 && j!=0)
	    {
	  //     if (i==N/2 && j==N/2)
// 		{
// 		  cout << ((*E1[pos])*(*E1[npos])).trace()  << endl;
// 		  cout << ((*E2[pos])*(*E2[npos])).trace()  << endl;
// 		  cout << ((*pi[pos])*(*pi[npos])).trace()  << endl << endl;
// 		  cout << ((*A1[pos])*(*A1[npos])).trace()  << endl;
// 		  cout << ((*A2[pos])*(*A2[npos])).trace()  << endl;
// 		  cout << ((*phi[pos])*(*phi[npos])).trace()  << endl << endl;
// 		}

	      if(omega2!=0)
		{
		  nkt1 = 1./sqrt(omega2)/static_cast<double>(N*N) * ( 1./((it-0.5)*dtau)*( (((*E1[pos])*(*E1[npos])).trace()).real() 
											   + (((*E2[pos])*(*E2[npos])).trace()).real() )
 		 					     );

		  nkt2 = 1./sqrt(omega2)/static_cast<double>(N*N) * (((it-0.5)*dtau) * ( (((*pi[pos])*(*pi[npos])).trace()).real() ));
		  
		  
		  nkt3 = sqrt(omega2)/static_cast<double>(N*N) * ( (it)*dtau * ( (((*A1[pos])*(*A1[npos])).trace()).real() 
									       + (((*A2[pos])*(*A2[npos])).trace()).real())
								   );
		  
		  nkt4 = sqrt(omega2)/static_cast<double>(N*N) * ( 1./((it)*dtau) * ( (((*phi[pos])*(*phi[npos])).trace()).real()));
		  
		  nkt5 = 1./static_cast<double>(N*N) * (complex<double>(0.,1.)*((*E1[pos])*(*A1[npos])-(*A1[pos])*(*E1[npos])+
										(*E2[pos])*(*A2[npos])-(*A2[pos])*(*E2[npos])).trace()).real();
		 
		  nkt6 = 1./static_cast<double>(N*N) * (complex<double>(0.,1.)*((*pi[pos])*(*phi[npos])-(*phi[pos])*(*pi[npos])).trace()).real();
	
		  
		  if(param->getRunWithkt()==1)
		    {
		      nkt1 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		      nkt2 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		      nkt3 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		      nkt4 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		      nkt5 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		      nkt6 *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
		    }
		}
	      else
		{
		  nkt1 = 0.;
		  nkt2 = 0.;
		  nkt3 = 0.;
		  nkt4 = 0.;
		  nkt5 = 0.;
		  nkt6 = 0.;
		}

	      dNdeta1 += nkt1;
	      dNdeta2 += nkt2;
	      dNdeta3 += nkt3;
	      dNdeta4 += nkt4;
	      dNdeta5 += nkt5;
	      dNdeta6 += nkt6;
	          	      
	      nkt = nkt1+nkt2+nkt3+nkt4+nkt5+nkt6; 
	      //  nktNoMixedTerms = nkt1+nkt2+nkt3+nkt4; 

	      dNdeta += nkt; // total multiplicity
	      //     dNdetaNoMixedTerms += nktNoMixedTerms; // total multiplicity
	      
	      nkxky[i][j] = nkt;
	      
	      // anglePhi = atan2(ky,kx); 
	      // if(anglePhi<0)
	      // 	anglePhi += 2.*Pi;

	      //	      cout << "kx=" << kx << ", ky=" << ky << ",phi=" << anglePhi << endl;


	      // for(int ik=0; ik<bins; ik++)
	      // 	{
	      // 	  if (abs(sqrt(kt2))>ik*dkt && abs(sqrt(kt2))<=(ik+1)*dkt)
	      // 	    {
	      // 	      nk[ik]+=nkt /2./Pi/2./Pi*N*N;
	      // 	      //    nkNoMixedTerms[ik]+=nktNoMixedTerms /2./Pi/2./Pi*N*N;
	      // 	      counterk[ik]+=1; // number of entries in nk[ik]
	      // 	      for(int iphi=0; iphi<phiBins; iphi++)
	      // 		{
	      // 		  if (anglePhi > deltaPhi*iphi && anglePhi <= (iphi+1)*deltaPhi)
	      // 		    {
	      // 		      //			      cout << deltaPhi*iphi << " < " << anglePhi << " < " << (iphi+1)*deltaPhi << endl;
	      // 		      //  if (counter[ik][iphi]>50)
	      // 		      //		continue;
	      // 		      n[ik][iphi]+=(nkt) /2./Pi/2./Pi*N*N;
	      // 		      n2[ik][iphi]+=nkt/dkt/2/Pi/sqrt(kt2);
	      // 		      //   nNoMixedTerms[ik][iphi]+=(nktNoMixedTerms) /2./Pi/2./Pi*N*N;
	      // 		      // change normalization to continuum one 
	      // 		      counter[ik][iphi]+=1; // number of entries in n[ik][iphi]
	      // 		    }
	      // 		}
	      // 	    }
	      // 	}
	    }
	}
    }
  
  int i, j;
  double latkx, latky;
  double fracX, fracY;
  for(int ik=0; ik<bins; ik++)
    {
      k = ik*dkt;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  anglePhi = deltaPhi*iphi;
	  
	  kx = k * cos(anglePhi);
	  ky = k * sin(anglePhi);

	  i = floor(((kx)/2/Pi+0.5)*N+1e-10);
	  j = floor(((ky)/2/Pi+0.5)*N+1e-10);

	  latkx = (2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(N)));
	  latky = (2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(N)));

	  //	  cout << "k=" << k << ", phi=" << anglePhi << ", kx=" << kx << ", ky=" << ky << ", lattice kx=" << latkx << ", lattice ky=" << latky << endl; 
	  
	  fracX = (kx - latkx)/(2*Pi/static_cast<double>(N));
	  fracY = (ky - latky)/(2*Pi/static_cast<double>(N));
	  
	  if(i+1<N && j+1 < N)
	    n[ik][iphi] = ((1.-fracX)*(1.-fracY)*nkxky[i][j]
			   + (fracX)*(1.-fracY)*nkxky[i+1][j]
			   + (1.-fracX)*(fracY)*nkxky[i][j+1]
			   + (fracX)*(fracY)*nkxky[i+1][j+1]) /2./Pi/2./Pi*N*N;//dkt/2/Pi/sqrt(kx*kx+ky*ky);
	  else
	    n[ik][iphi] = 0.;

	  if(k==0)
	    n[ik][iphi] = 0.;
	    
	  
	  //  cout << n[ik][iphi] << " " << nkxky[i][j] << " " << nkxky[i+1][j] << " " << nkxky[i][j+1] << " " << nkxky[i+1][j+1] << endl;


	}
    }

  double m,P;
  m=param->getJacobianm(); // in GeV
  P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
  double result, fullResult, fullResult2;
  fullResult=0.;
  fullResult2=0.;
  
   // for(int iphi=0; iphi<phiBins; iphi++)
   //   {
   //     cout << iphi*deltaPhi << " " << counter[1][iphi] << " " << counter[10][iphi] << " " << counter[20][iphi] << endl; 
   //   }
  
  for(int ik=1; ik<bins; ik++)
    {
      if(counterk[ik]>0)
	{
	  nk[ik] = nk[ik]/static_cast<double>(counterk[ik]);
	  //	  nkNoMixedTerms[ik] = nkNoMixedTerms[ik]/static_cast<double>(counterk[ik]);
	  //	  cout <<  ik*dkt/a*0.1973269718 << " " << static_cast<double>(counterk[ik]) << endl;
	}
      result = 0.;
      //  cout << "ik=" << ik << ", k =" <<  ik*dkt/a*0.1973269718<< endl;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  //  if(counter[ik][iphi]>0)
	  // {
	      //	      cout <<ik << " "<< iphi << " " << counter[ik][iphi] << " " << n[ik][iphi] << " " << n[ik][iphi]/static_cast<double>(counter[ik][iphi]) <<endl;
	  //  n[ik][iphi] = n[ik][iphi]/static_cast<double>(counter[ik][iphi]);
	      //	      nNoMixedTerms[ik][iphi] = nNoMixedTerms[ik][iphi]/static_cast<double>(counter[ik][iphi]);
	      result += n[ik][iphi] *deltaPhi;
	      //}
	}
      fullResult += result * (ik)*dkt*dkt;
      fullResult2 += nk[ik] * 2.*Pi *(ik+0.5)*dkt*dkt;
  //cout << "nk=" << nk[ik]*2.*Pi << ", integrated n=" << result << endl;
}
  cout << "N=" << dNdeta << ", k integrated N=" << fullResult2 << endl;
  cout << "N=" << dNdeta << ", k and phi integrated N=" << fullResult << endl;



  // output dN/d^2k
  ofstream foutPhiMult(PhiMult_name.c_str(),ios::app); 
  if (it==1)
    {
      foutPhiMult << "3" << " " << bins << " " << phiBins << endl; // 3 is the number of times we read out. modify if needed.
    }
  for(int ik=1; ik<bins; ik+=4)
    {
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  foutPhiMult << it*dtau*a << " " << ik*dkt/a*0.1973269718 << " "  << iphi*deltaPhi << " " 
		      << n[ik][iphi]*a/0.1973269718*a/0.1973269718 << endl; //<< " " << nNoMixedTerms[ik][iphi]*a/0.1973269718*a/0.1973269718 << endl;
	}
    }
  foutPhiMult.close();

 
  // compute hadrons using fragmentation function 


  int hbins = 40;
  int kcounter=0;
  double Nhad = 0.;
  double Ehad = 0.;
      double Nh[hbins+1][phiBins], Ng;
      
      for (int ih=0; ih<=hbins; ih++)
	{
	  for (int iphi=0; iphi<phiBins; iphi++)
	    {
	      Nh[ih][iphi]=0.;
	    }
	}
      double z,frac;
      double mypt, kt;
      int ik;
      int steps=200;
      double dz = 0.95/static_cast<double>(steps);
      
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  for(int ih=0; ih<=hbins; ih++)
	    {
	      mypt = ih * dkt/a*0.1973269718; //(10./static_cast<double>(hbins));  // the hadron's p_T
	      
	      //    int countCont=0;
	      for (int iz=0; iz<steps; iz++)
		{
		  z = 0.05 + iz * dz;
		  
		  kt = mypt/z; // the gluon's k_T
		  
		  ik = static_cast<int>(floor(kt*a/0.1973269718/dkt-0.5+0.00000001));
		  
		  frac = (kt - (ik+0.5)*dkt/a*0.1973269718)/(dkt/a*0.1973269718);
		  
		  if(ik+1<bins && ik >=0)
		    {
		      Ng = ((1.-frac)*n[ik][iphi]+frac*n[ik+1][iphi])*a/0.1973269718*a/0.1973269718; // to make dN/d^2k_T fo k_T in GeV 
		      if(kt>2)
		      	Ng*=exp(-(kt-2)*0.5);
		      //      countCont++;
		    }
		  else 
		    Ng = 0.;
		  
		  //integrate over z
		  if(param->getUsePseudoRapidity()==0)
		    {
		      if (z == 0.05 || z == 1.)
			{
			  Nh[ih][iphi] += 1./(z*z) * Ng * frag->kkp(7,1,z,kt) *dz *0.5; 
			}
		      else 
			{
			  Nh[ih][iphi] += 1./(z*z) * Ng * frag->kkp(7,1,z,kt)*dz; 
			}
		    }
		  else
		    {
		      if (z == 0.05 || z == 1.)
			{
			  Nh[ih][iphi] += 1./(z*z) * Ng * 
			    2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
				  +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
				  +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))))
			    *dz *0.5; 
			}
		      else 
			{
			  Nh[ih][iphi] += 1./(z*z) * Ng * 
			    2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
				  +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
				  +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))))
			    *dz;
			}
		    }
		}
	      //	      cout << "Using " << countCont << " kt's for p_T= "<< mypt << endl;
	    }
	}
      // output dN/d^2k
      ofstream foutPhiMultHad(PhiMultHad_name.c_str(),ios::app); 
      if (it==1)
	{
	  foutPhiMultHad << "3" << " " << hbins << " " << phiBins << endl; // 3 is the number of times we read out. modify if needed.
	}
      for(int ih=0; ih<hbins; ih++)
	{
	  for(int iphi=0; iphi<phiBins; iphi++)
	    {
	      foutPhiMultHad << it*dtau*a << " " << ih * dkt/a*0.1973269718 << " "  << iphi*deltaPhi << " " 
			     << Nh[ih][iphi] << endl;
	    }
	}
      foutPhiMultHad.close();


  ofstream foutCorr(Corr_name.c_str(),ios::app); 
  foutCorr << it*dtau*a << " " << dNdeta1 << " " << dNdeta2 << " " << dNdeta3 << " " << dNdeta4 << " " << dNdeta5 << " " << dNdeta6 <<  " " << dNdeta  << " " << dNdetaNoMixedTerms << endl;
  foutCorr.close();
  
  
  for(int i=0; i<N*N; i++)
    {
      delete E1[i];
      delete E2[i];
      delete pi[i];
      delete A1[i];
      delete A2[i];
      delete phi[i];
    }
  
  delete[] E1;
  delete[] E2;
  delete[] pi;
  delete[] A1;
  delete[] A2;
  delete[] phi;

  delete gaugefix;
  cout << " done." << endl;
  param->setSuccess(1);
  return 1;

}

int Evolution::correlationsColor(Lattice *lat, Group *group, Parameters *param, int it)
{
  int N = param->getSize();
  int Nc = param->getNc();
  int npos, pos, posNew;
  double L = param->getL();
  double a = L/N; // lattice spacing in fm
  double Pi, kx, ky, kt2;
  double g = param->getg();
  Pi = param->getPi();
  int nn[2];
  nn[0] = N;
  nn[1] = N;
  double dtau = param->getdtau();
  double nkt, nkt1, nkt2, nkt3, nkt4, nkt5, nkt6, nktNoMixedTerms;
  int bins = 20;
  int binsmult = 100;
  double nmult[binsmult];
  int phiBins = 16;
  double n[bins][phiBins]; // |k_T|, phi array
  double n2p[bins][phiBins][bins][phiBins]; // |k_T|, phi, |p_T|, phip array
  double n2[bins][phiBins]; // |k_T|, phi array
  double nkxky[N][N]; // kx, ky array
  
  vector< vector <double> > *nkxkypxpy;
  nkxkypxpy = new vector< vector <double> >; // kx, ky, px, py array
  vector< double > *npxpy;
  npxpy = new vector< double >; // px, py array
    

  double nk[bins]; //|k_T| array
  double nNoMixedTerms[bins][phiBins]; //|k_T|, phi array
  double nkNoMixedTerms[bins]; //|k_T| array
  int counter[bins][phiBins]; 
  int counterk[bins];
  int countermult[binsmult];
  double dkt = 1./static_cast<double>(bins);//sqrt(8.)/static_cast<double>(bins);
  double dktmult = sqrt(8.)/static_cast<double>(binsmult);
  double dNdeta =0.;
  double dNdetaNoMixedTerms =0.;
  double dNdeta1 =0.;
  double dNdeta2 =0.;
  double dNdeta3 =0.;
  double dNdeta4 =0.;
  double dNdeta5 =0.;
  double dNdeta6 =0.;
  double anglePhi;
  double k;  
  double deltaPhi = 2.*Pi/static_cast<double>(phiBins);

  stringstream strmult_name;
  strmult_name << "multiplicityCorr" << param->getMPIRank() << ".dat";
  string mult_name;
  mult_name = strmult_name.str();

  stringstream strmult2_name;
  strmult2_name << "multiplicityCorrFromPhi" << param->getMPIRank() << ".dat";
  string mult2_name;
  mult2_name = strmult2_name.str();

  stringstream strmult3_name;
  strmult3_name << "multiplicityCorrFromPhiHadrons" << param->getMPIRank() << ".dat";
  string mult3_name;
  mult3_name = strmult3_name.str();

  stringstream strPhiMult_name;
  strPhiMult_name << "MultPhi" << param->getMPIRank() << ".dat";
  string PhiMult_name;
  PhiMult_name = strPhiMult_name.str();

  stringstream strPhi2ParticleMult_name;
  strPhi2ParticleMult_name << "MultPhi2Particle" << param->getMPIRank() << ".dat";
  string Phi2ParticleMult_name;
  Phi2ParticleMult_name = strPhi2ParticleMult_name.str();

  stringstream strPhiMultHad_name;
  strPhiMultHad_name << "MultPhiHadrons" << param->getMPIRank() << ".dat";
  string PhiMultHad_name;
  PhiMultHad_name = strPhiMultHad_name.str();
  
  stringstream strPhi2ParticleMultHad_name;
  strPhi2ParticleMultHad_name << "MultPhiHadrons2Particle" << param->getMPIRank() << ".dat";
  string Phi2ParticleMultHad_name;
  Phi2ParticleMultHad_name = strPhi2ParticleMultHad_name.str();

  for(int ik=0; ik<binsmult; ik++)
    {
      nmult[ik] = 0.;
      countermult[ik]=0;
    }
 

  cout << "Measuring correlations... " << endl;


  // fix transverse Coulomb gauge
  GaugeFix *gaugefix;
  gaugefix = new GaugeFix(nn);
  
  double maxtime;
  if ( param->getInverseQsForMaxTime() == 1 )
    {
      maxtime = 1./param->getAverageQs()*hbarc;
      cout << "maximal evolution time = " << maxtime << " fm" << endl; 
    }
  else
    {
      maxtime = param->getMaxtime(); // maxtime is in fm
    }

  int itmax = static_cast<int>(floor(maxtime/(a*dtau)+1e-10));
  gaugefix->FFTChi(lat,group,param,4000);



  // gauge is fixed
  Matrix U1(Nc,1.);
  Matrix U2(Nc,1.);
  Matrix U1dag(Nc,1.);
  Matrix U2dag(Nc,1.);
 
  Matrix **A1;
  A1 = new Matrix*[N*N];
  Matrix **A2;
  A2 = new Matrix*[N*N];
  Matrix **phi;
  phi = new Matrix*[N*N];
  
  Matrix **E1;
  Matrix **E2;
  Matrix **pi;
  E1 = new Matrix*[N*N];
  E2 = new Matrix*[N*N];
  pi = new Matrix*[N*N];
  
  Matrix **t1;
  t1 = new Matrix*[N*N];
  Matrix **t2;
  t2 = new Matrix*[N*N];
  Matrix **t3;
  t3 = new Matrix*[N*N];
  Matrix **t4;
  t4 = new Matrix*[N*N];
  Matrix **t5;
  t5 = new Matrix*[N*N];
  Matrix **t6;
  t6 = new Matrix*[N*N];
  Matrix **t7;
  t7 = new Matrix*[N*N];
  Matrix **t8;
  t8 = new Matrix*[N*N];

  Matrix zero(Nc,0.);

  for(int i=0; i<N*N; i++)
    {
      A1[i] = new Matrix(Nc,0.);
      A2[i] = new Matrix(Nc,0.);
      E1[i] = new Matrix(Nc,0.);
      E2[i] = new Matrix(Nc,0.);
      pi[i] = new Matrix(Nc,0.);
      phi[i] = new Matrix(Nc,0.);
      t1[i] = new Matrix(Nc,0.);
      t2[i] = new Matrix(Nc,0.);
      t3[i] = new Matrix(Nc,0.);
      t4[i] = new Matrix(Nc,0.);
      t5[i] = new Matrix(Nc,0.);
      t6[i] = new Matrix(Nc,0.);
      t7[i] = new Matrix(Nc,0.);
      t8[i] = new Matrix(Nc,0.);
    }
  
 //version that determines the exact log of U1 and U2:
 for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  U1 = lat->cells[pos]->getUx();
	  U2 = lat->cells[pos]->getUy();
	  
	  U1.logm();
	  U2.logm();

	  *A1[pos] = complex<double>(0.,-1.)*U1;
	  *A2[pos] = complex<double>(0.,-1.)*U2;
	}
    }

  double g2mu2A, g2mu2B, gfactor, alphas, Qs;
  double c = param->getc();
  double muZero = param->getMuZero();
  // determine running coupling value
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  
	  if(param->getRunningCoupling())
	    {
	  
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2A = lat->cells[pos]->getg2mu2A();
		}
	      else 
		g2mu2A = 0;
	      
	      if(pos>0 && pos<(N-1)*N+N-1)
		{
		  g2mu2B = lat->cells[pos]->getg2mu2B();
		}
	      else
		g2mu2B = 0;
	      
	      if(param->getRunWithQs()==2)
		{
		  if(g2mu2A > g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==0)
		{
		  if(g2mu2A < g2mu2B)
		    Qs = sqrt(g2mu2A*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		  else
		    Qs = sqrt(g2mu2B*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}
	      else if(param->getRunWithQs()==1)
		{
		  Qs = sqrt((g2mu2A+g2mu2B)/2.*param->getQsmuRatio()*param->getQsmuRatio()/a/a*0.1973269718*0.1973269718*param->getg()*param->getg());
		}

	      
	      if ( param->getRunWithLocalQs() == 1)
		{
		  // 3 flavors
		  alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*Qs/0.2,2./c),c)));
		  gfactor = g*g/(4.*PI*alphas);
		  //		  cout << "c Qs=" << Qs << endl;
		  // cout << alphas << endl;
		  // run with the local (in transverse plane) coupling
		}
	      else
		{
		  if ( param->getRunWithQs() == 0 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsmin()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 1 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQsAvg()/0.2,2./c),c)));
		  else if ( param->getRunWithQs() == 2 )
		    alphas = 4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*param->getAverageQs()/0.2,2./c),c)));
		  
		  gfactor = g*g/(4.*PI*alphas);
		}
	    }
	  else
	    gfactor = 1.;

	  if(param->getRunWithkt()==0)
	    {
	      //	      cout << "gfactor=" << gfactor << endl;
	      *E1[pos] = lat->cells[pos]->getE1()*sqrt(gfactor); // replace one of the 1/g in the lattice E^i by the running one
	      *E2[pos] = lat->cells[pos]->getE2()*sqrt(gfactor); // "
	      *pi[pos] = lat->cells[pos]->getpi()*sqrt(gfactor); // replace the only 1/g by the running one (physical pi goes like 1/g, like physical E^i)
	      *A1[pos] = (*A1[pos])*sqrt(gfactor); //"
	      *A2[pos] = (*A2[pos])*sqrt(gfactor); // "
	      *phi[pos] = lat->cells[pos]->getphi()*sqrt(gfactor); // replace the only 1/g by the running one (physical pi goes like 1/g, like physical E^i)
	      //	      cout << "c " << sqrt(gfactor) << endl;
	    }
	  else
	    {
	      *E1[pos] = lat->cells[pos]->getE1(); 
	      *E2[pos] = lat->cells[pos]->getE2(); 
	      *pi[pos] = lat->cells[pos]->getpi(); 
	      *phi[pos] = lat->cells[pos]->getphi(); 
	    }
	}
    }

  // do Fourier transforms
  fft->fftn(A1,A1,nn,2,1);
  fft->fftn(A2,A2,nn,2,1);
  fft->fftn(phi,phi,nn,2,1);

  fft->fftn(E1,E1,nn,2,1);
  fft->fftn(E2,E2,nn,2,1);
  fft->fftn(pi,pi,nn,2,1);


  for(int ik=0; ik<bins; ik++)
    {
      nk[ik] = 0.;
      nkNoMixedTerms[ik] = 0.;
      counterk[ik] = 0;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  n[ik][iphi] = 0.;
	  n2[ik][iphi] = 0.;
	  nNoMixedTerms[ik][iphi] = 0.;
	  counter[ik][iphi]=0;
	}
    }


  // leave out the first cell to make it symmetric 
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{
	  pos = i*N+j;
	  npos = (N-i)*N+(N-j);
	  
	  kx = 2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(N));
	  ky = 2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(N));
	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));
	  
	  // i=0 or j=0 have no negative k_T value available

	  if(i!=0 && j!=0)
	    {
	      if(kt2!=0)
		{
		  *t1[pos] = (1./((it-0.5)*dtau)*1./sqrt(kt2)) * ( (*E1[pos])*(*E1[npos]) + (*E2[pos])*(*E2[npos]) );
		  *t2[pos] = ((it*dtau)*sqrt(kt2)) * ( (*A1[pos])*(*A1[npos]) + (*A2[pos])*(*A2[npos]) );
		  *t3[pos] = (((it-0.5)*dtau)*1./sqrt(kt2)) * ( (*pi[pos])*(*pi[npos]) );
		  *t4[pos] = (sqrt(kt2)/((it)*dtau)) * ( (*phi[pos])*(*phi[npos]) );
		  *t5[pos] = ( complex<double>(0.,1.)  ) * ( (*E1[pos])*(*A1[npos]) + (*E2[pos])*(*A2[npos]) );
		  *t6[pos] = ( complex<double>(0.,-1.) ) * ( (*A1[pos])*(*E1[npos]) + (*A2[pos])*(*E2[npos]) );
		  *t7[pos] = ( complex<double>(0.,1.)  ) * ( (*pi[pos])*(*phi[npos]) );
		  *t8[pos] = ( complex<double>(0.,-1.) ) * ( (*phi[pos])*(*pi[npos]) );
 		}
	      else
		{
		  *t1[pos] = zero;
		  *t2[pos] = zero;
		  *t3[pos] = zero;
		  *t4[pos] = zero;
		  *t5[pos] = zero;
		  *t6[pos] = zero;
		  *t7[pos] = zero;
		  *t8[pos] = zero;
		}
	    }
	}
    }

  int posk;
  int posp;
  double px, py, pt2;

  nkxkypxpy->clear();
  for(int ik=0; ik<N; ik++)
    {
      for(int jk=0; jk<N; jk++)
	{
	  posk = ik*N+jk;
	 
	  kx = 2.*Pi*(-0.5+static_cast<double>(ik)/static_cast<double>(N));
	  ky = 2.*Pi*(-0.5+static_cast<double>(jk)/static_cast<double>(N));
	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));
	
	  // used in multiplicity above:
	  // nkt = 2./sqrt(omega2)/static_cast<double>(N*N) * ( g*g/(it*dtau)*( (((*E1[pos])*(*E1[npos])).trace()).real() 
	  // 								     + (((*E2[pos])*(*E2[npos])).trace()).real() )
	  // 						     + (it*dtau) * ( (((*pi[pos])*(*pi[npos])).trace()).real() ) );
	  
	  nkxky[ik][jk] = 1./pow((N),2.) * ( (*t1[posk]) + (*t2[posk]) + (*t3[posk]) + (*t4[posk]) // ).trace().real();
	  				     + (*t5[posk]) + (*t6[posk]) + (*t7[posk]) + (*t8[posk])).trace().real();

	  //nkxky[ik][jk] = 2./pow((N),2.) * ( (*t1[posk]) + (*t3[posk])   ).trace().real();

	  if(param->getRunWithkt()==1)
	    {
	      nkxky[ik][jk] *= g*g/(4.*PI*4.*PI/(9.* log(pow(pow(muZero/0.2,2./c) + pow(param->getRunWithThisFactorTimesQs()*sqrt(kt2)*0.1973269718/a/0.2,2./c),c))));
	    }
	}
    }

  for(int ik=0; ik<N; ik++)
    {
      for(int jk=0; jk<N; jk++)
	{
	  kx = 2.*Pi*(-0.5+static_cast<double>(ik)/static_cast<double>(N));
	  ky = 2.*Pi*(-0.5+static_cast<double>(jk)/static_cast<double>(N));
	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.));

	  for(int i=0; i<binsmult; i++)
	    {
	      if (abs(sqrt(kt2))>i*dktmult && abs(sqrt(kt2))<=(i+1)*dktmult)
		{
		  nmult[i]+=nkxky[ik][jk]*N*N/Pi/Pi/2./2.;
		  countermult[i]+=1; // number of entries in n[ik]
		}
	    }
	}
    }

  ofstream foutmultiplicity(mult_name.c_str(),ios::app); 
	
  for(int ik=0; ik<binsmult; ik++)
    {
      if(countermult[ik]>0)
	{
	  nmult[ik] = nmult[ik]/static_cast<double>(countermult[ik]);
	}

      if(it == itmax)
	{
	  foutmultiplicity << it*dtau*a << " " << ik*dktmult/a*0.1973269718 << " " 
		       << nmult[ik]*a/0.1973269718*a/0.1973269718 << endl;
	}
    }
  
  foutmultiplicity.close();


  // now bin in phi

  int i, j;
  int i2, j2;
  double latkx, latky;
  double latpx, latpy;
  double fracX, fracY;
  double fracXp, fracYp;
  double p;
  double anglePhip;

  for(int ik=0; ik<bins; ik++)
    {
      k = ik*dkt;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  anglePhi = deltaPhi*iphi;
	  
	  kx = k * cos(anglePhi);
	  ky = k * sin(anglePhi);

	  i = floor(((kx)/2./Pi+0.5)*N+1e-10);
	  j = floor(((ky)/2./Pi+0.5)*N+1e-10);

	  latkx = (2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(N)));
	  latky = (2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(N)));

	  //	  cout << "k=" << k << ", phi=" << anglePhi << ", kx=" << kx << ", ky=" << ky << ", lattice kx=" << latkx << ", lattice ky=" << latky << endl; 
	  fracX = (kx - latkx)/(2*Pi/static_cast<double>(N));
	  fracY = (ky - latky)/(2*Pi/static_cast<double>(N));
	  
	  //check factors of N and 2pi

	  if(i+1<N && j+1 < N)
	    {
	      n[ik][iphi] = ((1.-fracX)*(1.-fracY)*nkxky[i][j]
			     + (fracX)*(1.-fracY)*nkxky[i+1][j]
			     + (1.-fracX)*(fracY)*nkxky[i][j+1]
			     + (fracX)*(fracY)*nkxky[i+1][j+1]) /k/2/Pi *N*N/2/Pi/2/Pi;
	    }
	  else
	    n[ik][iphi] = 0.;

	  if(k==0)
	    n[ik][iphi] = 0.;
	}
    }

  // for(int ik=0; ik<bins; ik++)
  //   {
  //     for(int iphi=0; iphi<phiBins; iphi++)
  // 	{
  // 	  if (counter[ik][iphi]>0)
  // 	    n[ik][iphi]/=static_cast<double>(counter[ik][iphi]);
  // 	}
  //   }

  double m,P;
  m=param->getJacobianm(); // in GeV
  P=0.13+0.32*pow(param->getRoots()/1000.,0.115); //in GeV
  double result, fullResult, fullResult2;
  fullResult=0.;
  fullResult2=0.;
  
  for(int ik=1; ik<bins; ik++)
    {
      if(counterk[ik]>0)
	{
	  nk[ik] = nk[ik]/static_cast<double>(counterk[ik]);
	}
      result = 0.;
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  result += n[ik][iphi] *deltaPhi;
	}
      fullResult += result * (ik)*dkt*dkt;
      fullResult2 += nk[ik] * 2.*Pi *(ik+0.5)*dkt*dkt;
    }
  cout << "N=" << dNdeta << ", k integrated N=" << fullResult2 << endl;
  cout << "N=" << dNdeta << ", k and phi integrated N=" << fullResult << endl;


  if(it == itmax)
    {
      ofstream foutmultiplicity2(mult2_name.c_str(),ios::app); 
      
      for(int ik=0; ik<bins; ik++)
	{
	  result = 0;
	  for(int iphi=0; iphi<phiBins; iphi++)
	    {
	      result += n[ik][iphi]*deltaPhi;
	    }
	  
	  foutmultiplicity2 << it*dtau*a << " " << ik*dkt/a*0.1973269718 << " " 
			    << result*a/0.1973269718*a/0.1973269718 << " " 
			    << n[ik][0]*a/0.1973269718*a/0.1973269718 << endl;
	}
      
      
      foutmultiplicity2.close();
    }


  // output dN/d^2k
  ofstream foutPhiMult(PhiMult_name.c_str(),ios::app); 
  if (it==1)
    {
      foutPhiMult << "3" << " " << (bins-1) << " " << phiBins << endl; // 3 is the number of times we read out. modify if needed.
    }
  for(int ik=1; ik<bins; ik+=1)
    {
      for(int iphi=0; iphi<phiBins; iphi++)
	{
	  foutPhiMult << it*dtau*a << " " << ik*dkt/a*0.1973269718 << " "  << iphi*deltaPhi << " " 
		      << n[ik][iphi]*a/0.1973269718*a/0.1973269718 << endl; //<< " " << nNoMixedTerms[ik][iphi]*a/0.1973269718*a/0.1973269718 << endl;
	}
    }
  foutPhiMult.close();




  // // output d^2N/d^2kd^2p
  // ofstream foutPhi2ParticleMult(Phi2ParticleMult_name.c_str(),ios::app); 
  // if (it==1)
  //   {
  //     foutPhi2ParticleMult << "3" << " " << (bins)/4 << " " << phiBins << endl; // 3 is the number of times we read out. modify if needed.
  //   }
  // for(int ik=1; ik<bins; ik+=4)
  //   {
  //     for(int iphi=0; iphi<phiBins; iphi++)
  // 	{
  // 	  for(int ip=1; ip<bins; ip+=4)
  // 	    {
  // 	      for(int iphip=0; iphip<phiBins; iphip++)
  // 		{
  // 		  foutPhi2ParticleMult << it*dtau*a << " " << ik*dkt/a*0.1973269718 << " "  << iphi*deltaPhi << " " 
  // 				       << " " << ip*dkt/a*0.1973269718 << " "  << iphip*deltaPhi << " " 
  // 				       << n2p[ik][iphi][ip][iphip]*a/0.1973269718*a/0.1973269718*a/0.1973269718*a/0.1973269718 << endl; 
  // 		}
  // 	    }
  // 	}
  //   }
  // foutPhi2ParticleMult.close();
 
  // compute hadrons using fragmentation function 
  int hbins = 20;
  int kcounter=0;
  double Nhad = 0.;
  double Ehad = 0.;
  double Nh[hbins+1][phiBins], Ng;
  
  for (int ih=0; ih<=hbins; ih++)
    {
      for (int iphi=0; iphi<phiBins; iphi++)
	{
	  Nh[ih][iphi]=0.;
	}
    }
    
  double z,z2,frac,frac2;
  double mypt2, mypt, kt;
  int ik, ik2;
  int steps=40;
  double dht = 0.25;
  double dz = 0.95/static_cast<double>(steps);
  
      for(int iphi=0; iphi<phiBins; iphi++)
  	{
  	  for(int ih=0; ih<=hbins; ih++)
  	    {
  	      mypt = ih * dht; //(10./static_cast<double>(hbins));  // the hadron's p_T
	      
  	      //    int countCont=0;
  	      for (int iz=0; iz<=steps; iz++)
  		{
  		  z = 0.05 + iz * dz;
		  
  		  kt = mypt/z; // the gluon's k_T
		  
  		  // ik = static_cast<int>(floor(kt*a/0.1973269718/dkt-0.5+0.00000001));
		  
  		  // frac = (kt - (ik+0.5)*dkt/a*0.1973269718)/(dkt/a*0.1973269718);
	
		  ik = static_cast<int>(floor(kt*a/0.1973269718/dkt+0.00000001));
		  
  		  frac = (kt - (ik)*dkt/a*0.1973269718)/(dkt/a*0.1973269718);
		  
  		  if(ik+1<bins && ik >=0)
  		    {
  		      Ng = ((1.-frac)*n[ik][iphi]+frac*n[ik+1][iphi])*a/0.1973269718*a/0.1973269718; // to make dN/d^2k_T for k_T in GeV 
  		    }
  		  else 
  		    Ng = 0.;
		  
      		  //integrate over z
  		  if(param->getUsePseudoRapidity()==0)
  		    {
  		      if (z == 0.05 || z == 1.)
  			{
  			  Nh[ih][iphi] += 1./(z*z) * Ng * frag->kkp(7,1,z,kt) *dz *0.5; 
  			}
  		      else 
  			{
  			  Nh[ih][iphi] += 1./(z*z) * Ng * frag->kkp(7,1,z,kt)*dz; 
  			}
  		    }
  		  else
  		    {
  		      if (z == 0.05 || z == 1.)
  			{
  			  Nh[ih][iphi] += 1./(z*z) * Ng * 
  			    2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
  				  +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
  				  +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))))
  			    *dz *0.5; 
  			}
  		      else 
  			{
  			  Nh[ih][iphi] += 1./(z*z) * Ng * 
  			    2. * (frag->kkp(1,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.13957*0.13957/(mypt*mypt)))
  				  +frag->kkp(2,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.493667*0.493667/(mypt*mypt)))
  				  +frag->kkp(4,1,z,kt)*cosh(param->getRapidity())/(sqrt(pow(cosh(param->getRapidity()),2.)+0.938272*0.938272/(mypt*mypt))))
  			    *dz;
  			}
  		    }
  		}
	    }
	}

      // output dN/d^2k
      ofstream foutPhiMultHad(PhiMultHad_name.c_str(),ios::app); 
      if (it==1)
  	{
  	  foutPhiMultHad << "3" << " " << hbins << " " << phiBins << endl; // 3 is the number of times we read out. modify if needed.
  	}
      for(int ih=0; ih<hbins; ih++)
  	{
  	  for(int iphi=0; iphi<phiBins; iphi++)
  	    {
  	      foutPhiMultHad << it*dtau*a << " " << ih * dht << " "  << iphi*deltaPhi << " " 
  			     << Nh[ih][iphi] << endl;
  	    }
  	}
      foutPhiMultHad.close();


      ofstream foutdNdpt(mult3_name.c_str(),ios::out); 
      for (int ih=0; ih<=hbins; ih++)
	{
	  result = 0.;
	  for(int iphi=0; iphi<phiBins; iphi++)
	    {
	      result += Nh[ih][iphi]*deltaPhi;
	      //	      cout << result << " " << 2*M_PI* Nh[ih][iphi] << endl;
	    }
	  foutdNdpt << ih*dht << " " << result << endl;
	}
      foutdNdpt.close();


  
  for(int i=0; i<N*N; i++)
    {
      delete E1[i];
      delete E2[i];
      delete pi[i];
      delete A1[i];
      delete A2[i];
      delete phi[i];
    }
  
  delete[] E1;
  delete[] E2;
  delete[] pi;
  delete[] A1;
  delete[] A2;
  delete[] phi;

  delete gaugefix;
  cout << " done." << endl;
  param->setSuccess(1);
  return 1;

}

  
  
void Evolution::twoPointFunctionInK(Parameters *param, Lattice *lat, int ids)
{
  int size = param->getSize(); 
  int Nc = param->getNc();
  double Pi = param->getPi();
  int bins;
  bins = size/2;
  int nn[2];
  nn[0]=size;
  nn[1]=size;

  int nbin[1];
  nbin[0]=bins;

  Matrix UDag(Nc,0);

  Matrix ** U;
  U = new Matrix *[size*size];
  
  for (int i=0; i<size*size; i++)
    {
      U[i] = new Matrix(Nc,0);
      *U[i] = lat->cells[i]->getU();
    }
  
  double **C;
  C = new double *[bins]; 
  double count[bins]; // counts entries in bin i
 
  for (int i=0; i<bins; i++)
    {
      C[i] = new double;
      *C[i] = 0.;
      count[i]=0.;
    }

  double kRange = 2.*sqrt(2.);
  double step = kRange/static_cast<double>(bins);
  int position;

  fft->fftn(U,U,nn,2,1);

  int pos;
  double kt, kx, ky;
  for (int i=0; i<size*size; i++)
    {
      UDag = *U[i];
      UDag.conjg();
      *U[i] = UDag*(*U[i]);
    }

  double sum = 0.;

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  // continuum momentum from FFT
	  kx = 2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
	  ky = 2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
	  kt = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); // lattice momentum
	  
	  position = static_cast<int>(floor(kt/step));
	
	  if (position<bins)
	    {
	      *C[position] += real((*U[pos]).trace())*kt*kt;
	      count[position]++;
	    }


// 	  if (position<bins)
// 	    {
// 	      if(Nc==2)
// 		{
// 		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(3))*kt*kt;
// 		  count[position]++;
// 		}
// 	      else if (Nc==3)
// 		{
// 		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))*kt*kt;
// 		  count[position]++;
// 		}
// 	    }
	}
    }

//   if (ids==0)
//     {
//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      if (Nc==2)
// 		{
// 		  sum+=(U[pos]->getRe(0) + U[pos]->getRe(3))/size/size/size/size;
// 		}
// 	      else if (Nc==3)
// 		{
// 		  sum+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))/size/size/size/size;
// 		}
// 	    }
// 	}  
//       cout << " sum=" << sum << endl;
//     }

  for(int i=0; i<bins; i++)
    {
      *C[i]/=static_cast<double>(count[i])*size*size; // divide by N^2, because sum_k C_k = N^2* N_c (instead of (2\pi)^2 N_c in continuum
    }

  ofstream fout("k-corr.dat",ios::app); 
  for (int j=0; j<bins; j++)
    {
      fout << (j*step+step/2.) << " "  << *C[j] << endl; // not scaled by g^2 mu (was k->k/g2mu  and C[j]->C[j]/g2mu/g2mu
    }
  fout << endl;
  fout.close();

  for (int i=0; i<size*size; i++)
    {
      delete U[i];
    }
  
  delete [] C;
  delete [] U;
}




  

