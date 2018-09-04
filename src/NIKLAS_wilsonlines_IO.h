#ifndef NIKLAS_output_wilsonlines_h
#define NIKLAS_output_wilsonlines_h



namespace WILSONLINES_IO
{

    void GenerateFileName(INT format, std::string &fname1, std::string &fname2, INT MPI_ID, INT MPI_SIZE, DOUBLE yeff, INT yeff_ID)
    {
        std::string fname = IO::OutDirectory+"/Wilsonlines_yeff"+to_string(yeff_ID);
        
        // check if this ID already exists //
        INT notfund=1;
        INT i=0;
        
        while(notfund)
        {
            INT offset = 2*i*MPI_SIZE;
            INT ID1 = offset + 2*MPI_ID ;
            INT ID2 = offset + 2*MPI_ID+1;
            
            std::string temp_file1,temp_file2;
            
            if(format==0)
            {
                temp_file1 = fname+"_"+to_string(ID1)+".txt";
                temp_file2 = fname+"_"+to_string(ID2)+".txt";
            }
            else if(format==1)
            {
                temp_file1 = fname+"_"+to_string(ID1)+".bin";
                temp_file2 = fname+"_"+to_string(ID2)+".bin";
            }
            else
            {
                std::cerr << "ERROR [GenerateFileName], data format not recognized (TXT,BIN): " << IO::DataFormat << std::endl;
                exit(0);
            }
            
            
            std::ifstream TestFile1(temp_file1.c_str());
            std::ifstream TestFile2(temp_file2.c_str());
            
            if( TestFile1.is_open() )
            {
                if( !TestFile2.is_open() )
                {
                    std::cerr << "# ERROR the second file (2) of the pair of files does not exist" << std::endl;
                    exit(0);
                }
                std::cout << "# files " << temp_file1 << " and " << temp_file2 << " already exists" << std::endl;
                
            }
            else
            {
                if( TestFile2.is_open() )
                {
                    std::cerr << "# ERROR the second file (2) EXISTS already, but the first does not " << std::endl;
                    exit(0);
                }
                std::cout << "# files " << temp_file1 << " and " << temp_file2 <<  " will be created" << std::endl;
                
                fname1=temp_file1;
                fname2=temp_file2;
                notfund=0;
            }
            i++;
        }
    }
    

    void GenerateFileName(std::string &fname1, std::string &fname2, INT MPI_ID, INT MPI_SIZE, DOUBLE yeff, INT yeff_ID)
    {
        std::string fname = IO::OutDirectory+"/Wilsonlines_yeff"+to_string(yeff_ID);
        
        // check if this ID already exists //
        INT notfund=1;
        INT i=0;
        
        while(notfund)
        {
            INT offset = 2*i*MPI_SIZE;
            INT ID1 = offset + 2*MPI_ID ;
            INT ID2 = offset + 2*MPI_ID+1;
            
            std::string temp_file1,temp_file2;
            
            if(IO::DataFormat=="TXT")
            {
                temp_file1 = fname+"_"+to_string(ID1)+".txt";
                temp_file2 = fname+"_"+to_string(ID2)+".txt";
            }
            else if(IO::DataFormat=="BIN")
            {
                temp_file1 = fname+"_"+to_string(ID1)+".bin";
                temp_file2 = fname+"_"+to_string(ID2)+".bin";
            }
            else
            {
                std::cerr << "ERROR [GenerateFileName], data format not recognized (TXT,BIN): " << IO::DataFormat << std::endl;
                exit(0);
            }
            
            
            std::ifstream TestFile1(temp_file1.c_str());
            std::ifstream TestFile2(temp_file2.c_str());
            
            if( TestFile1.is_open() )
            {
                if( !TestFile2.is_open() )
                {
                    std::cerr << "# ERROR the second file (2) of the pair of files does not exist" << std::endl;
                    exit(0);
                }
                std::cout << "# files " << temp_file1 << " and " << temp_file2 << " already exists" << std::endl;
                
            }
            else
            {
                if( TestFile2.is_open() )
                {
                    std::cerr << "# ERROR the second file (2) EXISTS already, but the first does not " << std::endl;
                    exit(0);
                }
                std::cout << "# files " << temp_file1 << " and " << temp_file2 <<  " will be created" << std::endl;
                
                fname1=temp_file1;
                fname2=temp_file2;
                notfund=0;
            }
            i++;
        }
    }
    
    
    INT Indx(INT ix, INT iy, INT N)
    {
        return ix*N+iy;
    }
    
    void ToCoords(INT ix, INT iy, DOUBLE &x, DOUBLE &y, DOUBLE L, DOUBLE a)
    {
        x = -L/2.+a*ix;
        y = -L/2.+a*iy;
    }
    
    
    Lattice *ReadWilsonLines(Parameters *param, std::string fname, INT &N, INT &Nc, DOUBLE &L, DOUBLE &a, INT MPI_ID, DOUBLE &yeff)
    {
        std::ifstream InStream(fname.c_str());
        DOUBLE re,im;
        Lattice *temp_lat;

        if(InStream.is_open())
        {
            // GRID CONFIGURATION //
            std::string str;
            char dummy;
            std::getline(InStream, str);
            std::stringstream ss(str);
            
            ss >> dummy;
            ss >> N;
            ss >> Nc;
            ss >> L;
            ss >> a;
            ss >> yeff;
            
            std::cout << "# Size is " << N << ", Nc " << Nc << ", length is " << L << ", a is [fm] " << a << std::endl;
            
            param = new Parameters;
            param->setSize(N);
            param->setNc(Nc);
            param->setL(L);
            param->setMPIRank(MPI_ID);
            temp_lat=new Lattice(param, Nc, N);      // first argument is the color representation, the second is the number of lattice sites per dimension //
            
            Matrix UBuffer(Nc);
            
            for(INT i=0; i<N*N; i++)
            {
                std::getline(InStream, str);
                std::stringstream buffer(str);
                
                DOUBLE x,y;
                buffer >> x;
                buffer >> y;

                for(INT j=0; j<3; j++)
                {
                    for(INT k=0; k<3; k++)
                    {
                        buffer >> re;
                        buffer >> im;
                        
                        UBuffer.set(k,j,COMPLEX(re,im));
                    }
                }
                
                temp_lat->cells[i]->setU(UBuffer);
            }
        }
    
        InStream.close();
        
        return temp_lat;
    }

    void PrintWilsonLinesBinary(Lattice *lat, Parameters *param, std::string fname1, std::string fname2, DOUBLE yeff )
    {
        INT N = param->getSize();
        INT Nc = param->getNc();
        DOUBLE L = param->getL();
        DOUBLE a = L/DOUBLE(N);          // lattice spacing in fm
        
        std::cout << "# PRINTING WILSON LINES, yeff = " << yeff << std::endl;
        
        std::ofstream Outfile1, Outfile2;

        Outfile1.precision(10);
        Outfile2.precision(10);
        
        Outfile1.open(fname1.c_str(),  ios::out | ios::binary);
        Outfile2.open(fname2.c_str(),  ios::out | ios::binary);
        
        DOUBLE temp=yeff;
        
        // print header ------------- //
        Outfile1.write((char *) &N ,sizeof(int));
        Outfile1.write((char *) &Nc ,sizeof(int));
        Outfile1.write((char *) &L ,sizeof(double));
        Outfile1.write((char *) &a ,sizeof(double));
        Outfile1.write((char *) &temp ,sizeof(double));
  
        Outfile2.write((char *) &N ,sizeof(int));
        Outfile2.write((char *) &Nc ,sizeof(int));
        Outfile2.write((char *) &L ,sizeof(double));
        Outfile2.write((char *) &a ,sizeof(double));
        Outfile2.write((char *) &temp ,sizeof(double));
        
        DOUBLE *val1=new DOUBLE[2];
        DOUBLE *val2=new DOUBLE[2];
        
        for(INT ix=0; ix<N; ix++)
        {
            for(INT iy=0; iy<N; iy++)
            {
                for(INT a=0; a<3; a++)
                {
                    for(INT b=0; b<3; b++)
                    {
                        val1[0] = (lat->cells[Indx(ix,iy,N)]->getU()).getRe(a*Nc+b);
                        val1[1] = (lat->cells[Indx(ix,iy,N)]->getU()).getIm(a*Nc+b);
                        val2[0] = (lat->cells[Indx(ix,iy,N)]->getU2()).getRe(a*Nc+b);
                        val2[1] = (lat->cells[Indx(ix,iy,N)]->getU2()).getIm(a*Nc+b);
                        
                       // std::cout << COMPLEX(val1[0],val1[1]) << std::endl;
                        
                        Outfile1.write((char *) val1 ,2*sizeof(DOUBLE));
                        Outfile2.write((char *) val2 ,2*sizeof(DOUBLE));
                    }
                }
            }
        }
        
        if(Outfile1.good()==false || Outfile2.good()==false)
        {
            std::cerr << "#CRTICAL ERROR -- BINARY OUTPUT OF VECTOR CURRENTS FAILED" << std::endl;
            exit(1);
        }
        
        
        delete [] val1;
        delete [] val2;
        
        Outfile1.close();
        Outfile2.close();
    }
    
    
    void PrintWilsonLines(Lattice *lat, Parameters *param, std::string fname1, std::string fname2, DOUBLE yeff )
    {
        INT N = param->getSize();
        INT Nc = param->getNc();
        DOUBLE L = param->getL();
        DOUBLE a = L/DOUBLE(N);          // lattice spacing in fm
        
        std::cout << "# PRINTING WILSON LINES" << std::endl;
        
        std::ofstream Outfile1, Outfile2;
        Outfile1.open(fname1.c_str());
        Outfile2.open(fname2.c_str());
        
        // print header ------------- //
        Outfile1 << "#\t" << N << "\t" << Nc << "\t" << L << "\t" << a << "\t" << yeff << std::endl;
        Outfile2 << "#\t" << N << "\t" << Nc << "\t" << L << "\t" << a << "\t" << yeff << std::endl;
        
        DOUBLE x,y;
        
        for(INT ix=0; ix<N; ix++)
        {
            for(INT iy=0; iy<N; iy++)
            {
               // ToCoords(ix, iy, x, y, L, a);
                Outfile1 << (lat->cells[Indx(ix,iy,N)]->getU()).MatrixToString() << std::endl;
                Outfile2 << (lat->cells[Indx(ix,iy,N)]->getU2()).MatrixToString() << std::endl;
            }
        }
        
        Outfile1.close();
        Outfile2.close();
    }
    
    void CompareLattice(Lattice *lat1, Lattice *lat2, Parameters *param)
    {
        std::cout << "# COMPARING LATTICES" << std::endl;
        
        INT N = param->getSize();
        INT Nc = param->getNc();
        
        std::cout << "# SIZE should be " << N << ", Nc " << Nc << std::endl;
        
        for(INT ix=0; ix<N; ix++)
        {
            for(INT iy=0; iy<N; iy++)
            {
                for(INT a=0; a<3; a++)
                {
                    for(INT b=0; b<3; b++)
                    {
                        DOUBLE val1re = (lat1->cells[Indx(ix,iy,N)]->getU()).getRe(a*Nc+b);
                        DOUBLE val1im = (lat1->cells[Indx(ix,iy,N)]->getU()).getIm(a*Nc+b);
                        DOUBLE val2re = (lat2->cells[Indx(ix,iy,N)]->getU()).getRe(a*Nc+b);
                        DOUBLE val2im = (lat2->cells[Indx(ix,iy,N)]->getU()).getIm(a*Nc+b);
                        
                        
                        if( abs(val1re -val2re) > 1e-4 )
                        {
                            std::cerr << "# ERROR BOTH LATTICES DO NOT AGREE, Re at (ix,iy) = " << ix << " " << iy << ", (a,b) = " << a << " " << b << std::endl;
                            exit(0);
                        }
                        if( abs(val1im -val2im) > 1e-4 )
                        {
                            std::cerr << "# ERROR BOTH LATTICES DO NOT AGREE, Im at (ix,iy) = " << ix << " " << iy << ", (a,b) = " << a << " " << b << std::endl;
                            exit(0);
                        }
                        
                    }
                }
            }
        }
        
        std::cout << "# LATTICES ARE IDENTICAL" << std::endl;
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
