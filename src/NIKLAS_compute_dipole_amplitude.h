#ifndef NIKLAS_compute_dipole_amplitude_h
#define NIKLAS_compute_dipole_amplitude_h



namespace DIPOLE_AMPLITUDE
{
    INT N;
    INT Nc;
    DOUBLE L;
    DOUBLE a;
    
    void Init(Parameters *param)
    {
        std::cout << "# INTIALIZING DIPOLE AMPLITUDE" << std::endl;
        N = param->getSize();
        Nc = param->getNc();
        L = param->getL();
        a = L/DOUBLE(N);    // lattice spacing in fm
    }
    
    INT Indx(INT ix, INT iy)
    {
        return ix*N+iy;
    }
    
    void ToCoords(INT ix, INT iy, DOUBLE &x, DOUBLE &y)
    {
        x = -L/2.+a*ix;
        y = -L/2.+a*iy;
    }
    
    void Compute(DATA::GRID *DipoleAmplitude, Lattice *lat)
    {
        std::cout << "# COMPUTING DIPOLE AMPLITUDES FROM WILSON LINES" << std::endl;

        Matrix U_dagger(Nc);
        Matrix U(Nc);
        DOUBLE x1, y1, x2, y2;

        for(INT i=0; i<DipoleAmplitude->dataSize; i++)
        {
            INT *r=new INT[DipoleAmplitude->NCoords];
            DipoleAmplitude->Coords(i, r);

            // U^dagger(x1,y1)
            U_dagger = lat->cells[Indx(r[0],r[1])]->getU();
            U_dagger.conjg();

            // U(x2,y2)
            U = lat->cells[Indx(r[2],r[3])]->getU();

            // N = 1 - (1/Nc)Tr ( U^dagger U ) //
            COMPLEX val = 1.0 - ( (U_dagger*U).trace() / DOUBLE(Nc) );
            
            ToCoords(r[0],r[1],x1,y1);
            ToCoords(r[2],r[3],x2,y2);

            DipoleAmplitude->p[0][r[0]]=x1;
            DipoleAmplitude->p[1][r[1]]=y1;
            DipoleAmplitude->p[2][r[2]]=x2;
            DipoleAmplitude->p[3][r[3]]=y2;

            DipoleAmplitude->Set(r,val);

            delete [] r;
        }

        std::cout << "# DIPOLE AMPILTUDE COMPUTED " << std::endl;
    }
    
    
    
    
    void PrintDipoleAmplitude(DATA::GRID *DipoleAmplitude, std::string fname)
    {
        // Fix b, print as function of r //
        std::ofstream Outfile;
        Outfile.open(fname.c_str());
        
        DOUBLE *x=new DOUBLE[DipoleAmplitude->NCoords];
        INT *r=new INT[DipoleAmplitude->NCoords];
        DOUBLE bx,by,rx,ry;
        DOUBLE abs_b, abs_r;
        
        for(INT i=0; i<DipoleAmplitude->dataSize; i++)
        {
            // COMPUTE COORDINATES //
            DipoleAmplitude->Coords(i, r);
            for(INT j=0; j<DipoleAmplitude->NCoords; j++)
            {
                x[j]=DipoleAmplitude->p[j][r[j]];
            }
            
            // CONVERT INTO ABS AND RELATIVE COORDS //
            KINEMATICS::SGL_to_CTR(x[0],x[1],x[2],x[3],bx,by,rx,ry);
            
            abs_b=std::sqrt(bx*bx+by*by);
            abs_r=std::sqrt(rx*rx+ry*ry);
            
            // PLOT FOR r=0 //
            if(abs_b == 0.0)
            {
                COMPLEX val=DipoleAmplitude->Get(r);
                Outfile << abs_r << " " << std::real(val) << " " << std::imag(val) << std::endl;
            }
        }
        
        Outfile.close();
        delete [] r;
        delete [] x;
    }

    
    
    void PrintTest(Lattice *lat, std::string fname)
    {
        std::cout << "# PRINTING TEST CASE" << std::endl;
        
        std::ofstream Outfile;
        Outfile.open(fname.c_str());
        
        Matrix U_dagger(Nc);
        Matrix U(Nc);
        DOUBLE x,y;
        
        U_dagger = lat->cells[Indx(0,0)]->getU();
        U_dagger.conjg();
        
        for(INT ix=0; ix<N; ix++)
        {
            for(INT iy=0; iy<N; iy++)
            {
                U = lat->cells[Indx(ix,iy)]->getU();
                DOUBLE val = std::real( (U_dagger*U).trace()) / DOUBLE(Nc);
                
                ToCoords(ix, iy, x, y);
                
                Outfile << x << " " << y << " " << val << std::endl;
            }
            Outfile << std::endl;
        }
        
        Outfile.close();
    }
    
    
    
    
    

}






    
    
    
//    void Compute(DATA::GRID *DipoleAmplitude, Lattice *lat)
//    {
//        Matrix Udag(Nc);
//
//
//
//        for(INT ix=0; ix<N; ix++)
//        {
//            for(INT iy=0; iy<N; iy++)
//            {
//
//            }
//        }
//
//
    
        
        // output correlator
        //  Udag = lat->cells[N*N/2+N/2]->getU();
        //Udag.conjg();
        
        // stringstream strVdagV_name;
        // strVdagV_name << "VdagV" << param->getMPIRank() << ".dat";
        // string VdagV_name;
        // VdagV_name = strVdagV_name.str();
        
        // ofstream fout(VdagV_name.c_str(),ios::out);
        
        // for (int i=0; i<nn[0]; i++)
        //     {
        //       for (int j=0; j<nn[1]; j++)
        //         {
        //           pos = i*N+j;
        //           x = -L/2.+a*i;
        //           y = -L/2.+a*j;
        
        //           corr = real((Udag*lat->cells[pos]->getU()).trace())/static_cast<double>(Nc);
        
        //           //            fout << x << ", " << y << ", " << corr << ", " << endl;
        //         }
        //       //        fout << endl;
        //     }
        
        // fout.close();
        
        
        // // output U
        // stringstream strVOne_name;
        // strVOne_name << "V1-" << param->getMPIRank() << ".txt";
        // string VOne_name;
        // VOne_name = strVOne_name.str();
        
        // ofstream foutU(VOne_name.c_str(),ios::out);
        // foutU.precision(15);
        
        // for(int ix=0; ix<N; ix++)
        //   {
        //     for(int iy=0; iy<N; iy++) // loop over all positions
        //     {
        //       pos = ix*N+iy;
        //       foutU << ix << " " << iy << " "  << (lat->cells[pos]->getU()).MatrixToString() << endl;
        //     }
        //     foutU << endl;
        //   }
        // foutU.close();
        
        // cout<<"wrote " << strVOne_name.str() <<endl;
        
        // stringstream strVTwo_name;
        // strVTwo_name << "V2-" << param->getMPIRank() << ".txt";
        // string VTwo_name;
        // VTwo_name = strVTwo_name.str();
        
        // ofstream foutU2(VTwo_name.c_str(),ios::out);
        // foutU2.precision(15);
        // for(int ix=0; ix<N; ix++)
        //   {
        //     for(int iy=0; iy<N; iy++) // loop over all positions
        //     {
        //       pos = ix*N+iy;
        //       foutU2 << ix << " " << iy << " "  << (lat->cells[pos]->getU2()).MatrixToString() << endl;
        //     }
        //     foutU2 << endl;
        //   }
        // foutU2.close();
        
        // cout<<"wrote " << strVTwo_name.str() <<endl;
        
   
    
    
    
    
    
    
    
    
    
    
    
    
































































#endif
