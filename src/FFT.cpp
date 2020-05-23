// FFT.cpp is part of the IP-Glasma solver.
// Copyright (C) 2012 Bjoern Schenke.
// This version uses FFTW
#include "FFT.h"

//**************************************************************************
// FFT class.

//**************************************************************************

void FFT::fftnVector(vector<complex<double> > **data, vector<complex<double> > **outdata, const int nn[], const int isign)
{
    unsigned ntot = nn[0]*nn[1];
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

void FFT::fftnArray(complex<double> **data, complex<double>  **outdata, const int nn[], const int isign, const int mDim)
{
    unsigned ntot = nn[0]*nn[1];
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
void FFT::fftn(T **data, T **outdata, const int nn[], const int isign)
{
    unsigned ntot = nn[0]*nn[1];
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
void FFT::fftnMany(T **data, T **outdata, const int nn[], const int isign)
{
    unsigned ntot = nn[0]*nn[1];
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

}


void FFT::fftnComplex(complex<double> *data, complex<double> *outdata, const int nn[], const int isign)
{
    unsigned ntot = nn[0]*nn[1];
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



// Define specializations of the template:
template void FFT::fftn(Matrix **data,Matrix **outdata, const int nn[], const int isign);
template void FFT::fftnMany(Matrix **data,Matrix **outdata, const int nn[], const int isign);
