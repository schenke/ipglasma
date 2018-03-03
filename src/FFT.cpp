// FFT.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
// This version uses FFTW
#include "FFT.h"

//**************************************************************************
// FFT class.

//**************************************************************************

void FFT::fftnVector(vector<complex<double> > **data, vector<complex<double> > **outdata, const int nn[], const int ndim, const int isign)
{
    unsigned i1, i2rev, i3rev, ibit;
    unsigned ip2, ifp1, ifp2, k2, n;
    unsigned nprev = 1, nrem, ntot = nn[0]*nn[1];
    register unsigned i2, i3;
    int pos, newpos;
    vector<complex<double> >::iterator position;	
		
    int mDim = data[0]->size();
    // mDim is the size of the vector (how many rows)

    // ndim is the dimension of the FFT (always 2 here)

    //    cout << "mDim=" << mDim << endl;
    // for each component of the vector fill the input array for the FFT (resort as you fill in)
    for (int k=0; k<mDim; k++)
      {
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = real(data[pos]->at(k));
		input[newpos][1] = imag(data[pos]->at(k));
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		input[newpos][0] = real(data[pos]->at(k));
		input[newpos][1] = imag(data[pos]->at(k));
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		input[newpos][0] = real(data[pos]->at(k));
		input[newpos][1] = imag(data[pos]->at(k));
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = real(data[pos]->at(k));
		input[newpos][1] = imag(data[pos]->at(k));
	      }
	  }	
	
	if (isign==1)
	  fftw_execute(p);
	else
	  fftw_execute(pback);

	// if this is inverse transform, normalize.
	if(isign == -1 ) 
	  {
	    for(unsigned i=0; i<ntot; i++) 
	      {
		output[i][0]/=static_cast<double>(ntot);
		output[i][1]/=static_cast<double>(ntot);
	      }
	  }
	
	//	xo   ->  oo
	//      oo       ox	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		position = outdata[pos]->begin() + k;	
		*position = complex<double>(output[newpos][0],output[newpos][1]);
		//	cout << output[newpos][0] << " " << output[newpos][1] << endl;
		//cout << "data=" << data[pos]->at(k) << endl;
	      }
	  }
	//	oo   ->  xo
	//      ox       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		position = outdata[pos]->begin() + k;	
		*position = complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		position = outdata[pos]->begin() + k;	
		*position = complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	

	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		position = outdata[pos]->begin() + k;	
		*position = complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
      }
}

void FFT::fftnArray(complex<double> **data, complex<double>  **outdata, const int nn[], const int ndim, const int isign, const int mDim)
{
    unsigned i1, i2rev, i3rev, ibit;
    unsigned ip2, ifp1, ifp2, k2, n;
    unsigned nprev = 1, nrem, ntot = nn[0]*nn[1];
    register unsigned i2, i3;
    int pos, newpos;
		
    //    cout << "size=" << mDim << endl;
    // mDim is the size of the vector (how many rows)

    // ndim is the dimension of the FFT (always 2 here)

    //    cout << "mDim=" << mDim << endl;
    // for each component of the vector fill the input array for the FFT (resort as you fill in)
    for (int k=0; k<mDim; k++)
      {
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = real(data[pos][k]);
		input[newpos][1] = imag(data[pos][k]);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		input[newpos][0] = real(data[pos][k]);
		input[newpos][1] = imag(data[pos][k]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		input[newpos][0] = real(data[pos][k]);
		input[newpos][1] = imag(data[pos][k]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = real(data[pos][k]);
		input[newpos][1] = imag(data[pos][k]);
	      }
	  }	
	
	if (isign==1)
	  fftw_execute(p);
	else
	  fftw_execute(pback);

	// if this is inverse transform, normalize.
	if(isign == -1 ) 
	  {
	    for(unsigned i=0; i<ntot; i++) 
	      {
		output[i][0]/=static_cast<double>(ntot);
		output[i][1]/=static_cast<double>(ntot);
	      }
	  }
	
	//	xo   ->  oo
	//      oo       ox	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos][k]= complex<double>(output[newpos][0],output[newpos][1]);
		//	cout << output[newpos][0] << " " << output[newpos][1] << endl;
		//cout << "data=" << data[pos][k] << endl;
	      }
	  }
	//	oo   ->  xo
	//      ox       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		outdata[pos][k]= complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		outdata[pos][k]= complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	

	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos][k]= complex<double>(output[newpos][0],output[newpos][1]);
	      }
	  }
      }
}



// Performs Fast Fourier Transform of any object of class "T" (matrix or something else) using a wrapper for FFTW
// This routine takes data as a function of -x_max/2 to x_max/2 and returns it ordered similarly - no need to resort before or after!
template <class T>
void FFT::fftn(T **data, T **outdata, const int nn[], const int ndim, const int isign)
{
    unsigned i1, i2rev, i3rev, ibit;
    unsigned ip2, ifp1, ifp2, k2, n;
    unsigned nprev = 1, nrem, ntot = ntot = nn[0]*nn[1];
    register unsigned i2, i3;
    int mDim, pos, newpos; // matrix dimension
    mDim = data[0]->getNDim();
    
    // mDim is the size of the matrix (how many rows)

    mDim*=mDim;
    // ndim is the dimension of the FFT (always 2 here)
    
    // for each component of the matrix fill the input array for the FFT (resort as you fill in)

    //    cout << "mDim=" << mDim << endl;
    for (int k=0; k<mDim; k++)
      {
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }	
	
	if (isign==1)
	  fftw_execute(p);
	else
	  fftw_execute(pback);
	
	// if this is inverse transform, normalize.
	if(isign == -1 ) 
	  {
	    for(unsigned i=0; i<ntot; i++) 
	      {
		output[i][0]/=static_cast<double>(ntot);
		output[i][1]/=static_cast<double>(ntot);
	      }
	  }
	
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	

	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
      }
}


template <class T>
void FFT::fftnMany(T **data, T **outdata, const int nn[], const int ndim, const int isign)
{
    unsigned i1, i2rev, i3rev, ibit;
    unsigned ip2, ifp1, ifp2, k2, n;
    unsigned nprev = 1, nrem, ntot = ntot = nn[0]*nn[1];
    register unsigned i2, i3;
    int mDim, pos, newpos; // matrix dimension
    mDim = data[0]->getNDim();
    
    // mDim is the size of the matrix (how many rows)

    mDim*=mDim;
    // ndim is the dimension of the FFT (always 2 here)
    
    // for each component of the matrix fill the input array for the FFT (resort as you fill in)
    
    //    cout << "mDim=" << mDim << endl;

   for (int k=0; k<mDim; k++)
      {
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2+k*ntot;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j+k;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j+k*ntot;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2+k*ntot;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2+k*ntot;
		input[newpos][0] = data[pos]->getRe(k);
		input[newpos][1] = data[pos]->getIm(k);
	      }
	  }	
      }
   
    if (isign==1)
      fftw_execute(pmany);
    else
      fftw_execute(pmanyback);

    // if this is inverse transform, normalize.
    if(isign == -1 ) 
      {
	for(unsigned i=0; i<ntot*9; i++) 
	  {
	    output[i][0]/=static_cast<double>(ntot);
	    output[i][1]/=static_cast<double>(ntot);
	  }
      }
   
   for (int k=0; k<mDim; k++)
     {

	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2+k*ntot;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j+k*ntot;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2+k*ntot;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	

	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2+k*ntot;
		outdata[pos]->setRe(k,output[newpos][0]);
		outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
      }

   //------


    // for (int k=0; k<mDim; k++)
    //   {
    // 	for (int i=0; i<nn[0]; i++)
    // 	  {
    // 	    for (int j=0; j<nn[1]; j++)
    // 	      {
    // 		pos = i*nn[1]+j;
    // 		newpos = i*nn[1]+j+k*ntot;
    // 		input[newpos][0] = data[pos]->getRe(k);
    // 		input[newpos][1] = data[pos]->getIm(k);
    // 	      }
    // 	  }
    //   }
    // if (isign==1)
    //   fftw_execute(pmany);
    // else
    //   fftw_execute(pmanyback);
    
    // // if this is inverse transform, normalize.
    // if(isign == -1 ) 
    //   {
    // 	for(unsigned i=0; i<ntot*9; i++) 
    // 	  {
    // 	    output[i][0]/=static_cast<double>(ntot);
    // 	    output[i][1]/=static_cast<double>(ntot);
    // 	  }
    //   }
    
    // for (int k=0; k<mDim; k++)
    //   {
    // 	for (int i=0; i<nn[0]; i++)
    // 	  {
    // 	    for (int j=0; j<nn[1]; j++)
    // 	      {
    // 		pos = i*nn[1]+j;
    // 		newpos = i*nn[1]+j+k*ntot;
    // 		outdata[pos]->setRe(k,output[newpos][0]);
    // 		outdata[pos]->setIm(k,output[newpos][1]);
    // 	      }
    // 	  }
    //   }
}


void FFT::fftnComplex(complex<double> *data, complex<double> *outdata, const int nn[], const int ndim, const int isign)
{
    unsigned i1, i2rev, i3rev, ibit;
    unsigned ip2, ifp1, ifp2, k2, n;
    unsigned nprev = 1, nrem, ntot = ntot = nn[0]*nn[1];
    register unsigned i2, i3;
    int mDim, pos, newpos; // matrix dimension
    mDim = 1;
    
    // mDim is the size of the matrix (how many rows)

    mDim*=mDim;
    // ndim is the dimension of the FFT (always 2 here)
    
    // for each component of the matrix fill the input array for the FFT (resort as you fill in)

    //    cout << "mDim=" << mDim << endl;
    for (int k=0; k<mDim; k++)
      {
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = data[pos].real();
		input[newpos][1] = data[pos].imag();
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		input[newpos][0] = data[pos].real();
		input[newpos][1] = data[pos].imag();
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		input[newpos][0] = data[pos].real();
		input[newpos][1] = data[pos].imag();
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		input[newpos][0] = data[pos].real();
		input[newpos][1] = data[pos].imag();
	      }
	  }	
	
	if (isign==1)
	  fftw_execute(p);
	else
	  fftw_execute(pback);
	
	// if this is inverse transform, normalize.
	if(isign == -1 ) 
	  {
	    for(unsigned i=0; i<ntot; i++) 
	      {
		output[i][0]/=static_cast<double>(ntot);
		output[i][1]/=static_cast<double>(ntot);
	      }
	  }
	
	//	oo   ->  xo
	//      ox       oo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos]=complex<double>(output[newpos][0],output[newpos][1]);
		//	outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	xo   ->  oo
	//      oo       ox	
	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
		outdata[pos]=complex<double>(output[newpos][0],output[newpos][1]);
		//	outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	ox   ->  oo
	//      oo       xo	
	for (int i=nn[0]/2; i<nn[0]; i++)
	  {
	    for (int j=0; j<nn[1]/2; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
		outdata[pos]=complex<double>(output[newpos][0],output[newpos][1]);
		//	outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
	//	oo   ->  ox
	//      xo       oo	

	for (int i=0; i<nn[0]/2; i++)
	  {
	    for (int j=nn[1]/2; j<nn[1]; j++)
	      {
		pos = i*nn[1]+j;
		newpos = (i+nn[0]/2)*nn[1]+j-nn[1]/2;
		outdata[pos]=complex<double>(output[newpos][0],output[newpos][1]);
		//	outdata[pos]->setIm(k,output[newpos][1]);
	      }
	  }
      }
}


// void FFT::fftnComplex(complex<double> **data, const int nn[], const int ndim, const int isign)
// {
//     unsigned i1, i2rev, i3rev, ibit;
//     unsigned ip2, ifp1, ifp2, k2, n;
//     unsigned nprev = 1, nrem, ntot = nn[0]*nn[1];
//     register unsigned i2, i3;
//     int mDim=1;

//     // Initialize auxiliary object as 0.
//     complex<double> *temp = new complex<double>; 
//     *temp = 0.;
//     complex<double> *temp2 = new complex<double>;
//     *temp2 = 0.;
    
//     // temporary grid needed for resorting
//     complex<double>** dataTemp;
//     dataTemp = new complex<double>*[nn[0]*nn[1]/4];
    
//     int newpos, pos; 
//     for(int i=0; i<nn[0]*nn[1]/4; i++)
//       {
// 	dataTemp[i] = new complex<double>;
// 	*dataTemp[i] = 0.;
//       }	  
    
//     for (int i=0; i<nn[0]/2; i++)
//       {
// 	for (int j=0; j<nn[1]/2; j++)
// 	  {
// 	    pos = i*nn[1]+j;
// 	    newpos = i*nn[1]/2+j;
// 	    *dataTemp[newpos] = *data[pos];
// 	  }
//       }
//     for (int i=nn[0]/2; i<nn[0]; i++)
//       {
// 	for (int j=nn[1]/2; j<nn[1]; j++)
// 	  {
// 	    pos = i*nn[1]+j;
// 	    newpos = (i-nn[0]/2)*nn[1]+j-nn[1]/2;
// 	    *data[newpos] = *data[pos];
// 	  }
//       }
//     for (int i=0; i<nn[0]/2; i++)
//       {
// 	for (int j=0; j<nn[1]/2; j++)
// 	  {
// 	    pos = i*nn[1]/2+j;
// 	    newpos = (i+nn[0]/2)*nn[1]+nn[1]/2+j;
// 	    *data[newpos] = *dataTemp[pos];
// 	  }
//       }
//     for (int i=0; i<nn[0]/2; i++)
//       {
// 	for (int j=nn[1]/2; j<nn[1]; j++)
// 	  {
// 	    pos = i*nn[1]+j;
// 	    newpos = i*nn[1]/2+j-nn[1]/2;
// 	    *dataTemp[newpos] = *data[pos];
// 	  }
//       }
//     for (int i=nn[0]/2; i<nn[0]; i++)
//       {
// 	for (int j=0; j<nn[1]/2; j++)
// 	  {
// 	    pos = i*nn[1]+j;
// 	    newpos = (i-nn[0]/2)*nn[1]+j+nn[1]/2;
// 	    *data[newpos] = *data[pos];
// 	  }
//       }
//     for (int i=0; i<nn[0]/2; i++)
//       {
// 	for (int j=0; j<nn[1]/2; j++)
// 	  {
// 	    pos = i*nn[1]/2+j;
// 	    newpos = (i+nn[0]/2)*nn[1]+j;
// 	    *data[newpos] = *dataTemp[pos];
// 	  }
//       }
    
//     // Compute total number of complex values.
//     for (int idim = 0; idim < ndim; ++idim) ntot *= nn[idim];

//     for (int idim = ndim - 1; idim >= 0; --idim) 
//       {
// 	n = nn[idim];
// 	nrem = ntot / (n * nprev);
// 	ip2 = nprev * n;        // Unit step for next dimension.
// 	i2rev = 0;              // Bit reversed i2.

// 	// This is the bit reversal section of the routine.
// 	// Loop over current dimension.
// 	for (i2 = 0; i2 < ip2; i2 += nprev) 
// 	  {
// 	    if (i2 < i2rev)
// 	      // Loop over lower dimensions.
// 	      for (i1 = i2; i1 < i2 + nprev; ++i1)
// 		// Loop over higher dimensions.
// 		for (i3 = i1; i3 < ntot; i3 += ip2) 
// 		  {
// 		    i3rev = i3 + i2rev - i2;
// 		    *temp = *(data[i3]);
// 		    *(data[i3]) = *(data[i3rev]);
// 		    *(data[i3rev]) = *temp;
// 		  }
	    
// 	    ibit = ip2;
// 	    // Increment from high end of i2rev to low.
// 	    do {
// 	      ibit >>= 1;
// 	      i2rev ^= ibit;
// 	    } while (ibit >= nprev && !(ibit & i2rev));
// 	  }

// 	// Here begins the Danielson-Lanczos section of the routine.
// 	for (ifp1 = nprev; ifp1 < ip2; ifp1 <<= 1) 
// 	  {
// 	    ifp2 = ifp1 << 1;
// 	    //  Initialize for the trig. recurrence.
// 	    double alp = sin(M_PI/ (ifp2/nprev));
// 	    alp = 2*alp*alp;
// 	    double beta = -isign*sin(2*M_PI/(ifp2/nprev));
	    
// 	    complex<double> wp(alp,beta);
// 	    complex<double> w(1.0,0.0);
	    
// 	    // Loop by unit step in current dimension.
// 	    for (i3 = 0; i3 < ifp1; i3 += nprev) 
// 	      {
// 		// Loop over lower dimensions.
// 		for (i1 = i3; i1 < i3 + nprev; ++i1)
// 		  // Loop over higher dimensions.
// 		  for (i2 = i1; i2 < ntot; i2 += ifp2) 
// 		    {
// 		      // Danielson-Lanczos formula.
// 		      k2 = i2 + ifp1;
// 		      *temp  = *(data[k2]);
// 		      *temp2 = *(data[k2]);
// 		      *temp  *= w;
// 		      *temp2 *= -w;
// 		      *temp  += *(data[i2]);
// 		      *temp2 += *(data[i2]);
// 		      *(data[i2]) = *temp;
// 		      *(data[k2]) = *temp2;
// 		    }
// 		// Trigonometric recurrence
// 		w -= w*wp;
// 	      }
// 	  }
// 	nprev *= n;
//       }
//     // if this is inverse transform, normalize.
//     if(isign == -1 ) 
//       {
// 	for(unsigned i=0; i<ntot; i++) 
// 	  {
// 	    *(data[i]) /= ntot;
// 	  }
//       }



//     for(int i=0; i<nn[0]*nn[1]/4; i++)
//       {
// 	delete dataTemp[i];
//       }	

//     delete[] dataTemp;
//     delete temp;
//     delete temp2;
// }


// Define specializations of the template:
template void FFT::fftn(Matrix **data,Matrix **outdata, const int nn[], const int ndim, const int isign);
template void FFT::fftnMany(Matrix **data,Matrix **outdata, const int nn[], const int ndim, const int isign);
