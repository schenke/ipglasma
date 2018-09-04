//
//  integrated_crosssections.h
//  SubNucleon diffraction
//
//  Created by Niklas Mueller on 4/7/18.
//  Copyright © 2018 Heikki. All rights reserved.
//
#ifndef datagrids_h
#define datagrids_h
#define INT int
#define DOUBLE double
#define COMPLEX complex<double>



namespace DATA
{
    class GRID
    {
    public:
        GRID(INT NNCoords, INT *GridSize, string *IIdentifier)
        {
            // number of coordinates for grid //
            NCoords = NNCoords;
            
            // size of grid in each dimension //
            Size = new INT[NCoords];
            
            // coordinates in each dimension //
            p = new DOUBLE*[NCoords];
            
            // Identifier //
            Identifier = new string[NCoords];
            
            // total size of values //
            dataSize=1;
            
            for(INT i=0; i<NCoords; i++)
            {
                Size[i] = GridSize[i];
                Identifier[i] = IIdentifier[i];
                p[i] = new DOUBLE[GridSize[i]];
                dataSize *=GridSize[i];
            }
            
            F= new complex<DOUBLE>[dataSize];
        }
        
        ~GRID()
        {
            delete [] F;
        }
        
        void Assign_MPI_tasks(INT MPI_ID, INT MPI_threads)
        {
            MPI_lower = new INT[MPI_threads];
            MPI_upper = new INT[MPI_threads];
            
            INT single_MPI_range = INT(this->dataSize/MPI_threads);
            
            for(INT i=0;i<MPI_threads; i++)
            {
                if(i<MPI_threads-1)
                {
                    MPI_lower[i] = i*single_MPI_range;
                    MPI_upper[i] = (i+1)*single_MPI_range;
                }
                else
                {
                    MPI_lower[i] = i*single_MPI_range;
                    MPI_upper[i] = this->dataSize;
                }
            }
            
            if(MPI_ID==0)
            {
                for(INT i=0;i<MPI_threads; i++)
                {
                    std::cout << "# data range MPI " << i << " : " << MPI_lower[i] <<  " " << MPI_upper[i] << std::endl;
                }
            }
        }
        
        void MPI_clean()
        {
            delete [] MPI_lower;
            delete [] MPI_upper;
        }
        
        void Init(DOUBLE *lower, DOUBLE *upper)
        {
            // initializes an equidistant grid //
            for(INT i=0; i<NCoords; i++)
            {
                for(INT j=0; j<Size[i]; j++)
                {
                    p[i][j] = lower[i] + (DOUBLE(j)/DOUBLE(Size[i]-1)) * (upper[i]-lower[i]);
                }
            }
        }
        
        
        
        
        // takes the n-dimensional integral over all data points//
        complex<DOUBLE> IntegrateAll()
        {
            DOUBLE VolumeElement=1.0;
            for(INT i=0; i<NCoords; i++)
            {
                VolumeElement *= (this->p[i][1]-this->p[i][0]);
            }
            
            complex<DOUBLE> sum=0.;
            INT *ppi = new INT[NCoords];
            for(INT i=0; i<dataSize; i++)
            {
                Coords(i,ppi);
                complex<DOUBLE> val = this->Get(ppi);
                IntegralUpdate(ppi, val, sum, VolumeElement, Size, NCoords);
            }
            return sum;
        }
        
        
        GRID *CloneEmpty()
        {
            GRID *SubGrid = new GRID(this->NCoords, this->Size, this->Identifier);
            
            DOUBLE *lower=new DOUBLE[NCoords];
            DOUBLE *upper=new DOUBLE[NCoords];
            
            for(INT i=0; i<NCoords; i++)
            {
                lower[i]=this->p[i][0];
                upper[i]=this->p[i][Size[i]-1];
                std::cout << lower[i] << " " << upper[i] << std::endl;
            }
            
            
            SubGrid->Init(lower,upper);
            
            return SubGrid;
        }
        
        bool IsOnGrid(DOUBLE *x)
        {
            bool result=true;
            for(INT i=0; i<NCoords; i++)
            {
                if( x[i]<this->p[i][0] || x[i]>this->p[i][Size[i]-1] )
                {
                    result=false;
                }
            }
            return result;
        }
        
        
        
        
        
        GRID *CreateSubGRID_Integrate(INT IntegralDIM, INT *Mask)
        {
            INT NewDim = NCoords-IntegralDIM;          // new dimension of subgrid //
            INT *SSize = new INT[NewDim];
            INT *AntiSSize = new INT[IntegralDIM];
            string *IIdentifier = new string[NewDim];
            string *AntiIIdentifier = new string[IntegralDIM];
            
            // BUILDING SUBGRID //
            INT k=0,antik=0;
            for(INT i=0; i<NCoords; i++)
            {
                if(Mask[i]==1)
                {
                    SSize[k]=this->Size[i];
                    IIdentifier[k]=this->Identifier[i];
                    k++;
                }
                else if(Mask[i]==0)
                {
                    AntiSSize[antik]=this->Size[i];
                    AntiIIdentifier[antik]=this->Identifier[i];
                    antik++;
                }
            }
            
            GRID *SubGrid = new GRID(NewDim, SSize, IIdentifier);
            GRID *IntegralGrid = new GRID(IntegralDIM, AntiSSize, AntiIIdentifier);
            
            // Copy coordinates from parent grid //
            k=0; antik=0;
            for(INT i=0; i<NCoords; i++)
            {
                if(Mask[i]==1)
                {
                    if( SubGrid->Size[k] != this->Size[i])
                    {
                        std::cerr << "sizes of parent and sub grid do not match, subgrid: " << SubGrid->Size[k] << ", parent: " << this->Size[i] << std::endl;
                        exit(0);
                    }
                    for(INT j=0; j<SubGrid->Size[k]; j++)
                    {
                        SubGrid->p[k][j]=this->p[i][j];
                    }
                    k++;
                }
                else if(Mask[i]==0)
                {
                    if( IntegralGrid->Size[antik] != this->Size[i])
                    {
                        std::cerr << "sizes of parent and integral grid do not match, integral grid: " << IntegralGrid->Size[antik] << ", parent: " << this->Size[i] << std::endl;
                        exit(0);
                    }
                    for(INT j=0; j<IntegralGrid->Size[antik]; j++)
                    {
                        IntegralGrid->p[antik][j]=this->p[i][j];
                    }
                    antik++;
                }
            }
            
            
            // Compute Integral values //
            INT *r_sub = new INT[SubGrid->NCoords];
            INT *r_int = new INT[IntegralDIM];
            INT *r     = new INT[this->NCoords];
            
            INT IntegralSize=this->dataSize/SubGrid->dataSize;
            if(IntegralSize != IntegralGrid->dataSize)
            {
                std::cerr << "sizes of integral grid is " << IntegralGrid->dataSize << " vs. " << IntegralSize << std::endl;
                exit(0);
            }
            std::cout << "# integral size " << IntegralSize << " (total " << this->dataSize << ", subgrid " << SubGrid->dataSize << ")" << std::endl;
            
            for(INT i=0; i<SubGrid->dataSize; i++)
            {
                complex<DOUBLE> sum=0.0;
                SubGrid->Coords(i,r_sub);
                
                for(INT j=0; j<IntegralSize; j++)
                {
                    IntegralGrid->Coords(j,r_int);
                    CombineIndices(r_sub, r_int, r, Mask);
                    
                    complex<DOUBLE> val=this->Get(r);

                    DOUBLE weight=1.0;
                    for(INT m=0; m<IntegralDIM; m++)
                    {
                        // volume element //
                        weight *=(IntegralGrid->p[m][1]-IntegralGrid->p[m][0]);

                        // check if a boundary term //
                        if( r_int[m]==0 || r_int[m]==IntegralGrid->Size[m]-1 )
                        {
                            weight *=0.5;
                        }
                    }

                    sum += weight*val;
                }
                
                SubGrid->Set(r_sub,sum);
            }
            
            delete IntegralGrid;
            delete [] r_sub;
            delete [] r_int;
            delete [] r;
            delete [] SSize;
            delete [] AntiSSize;
            delete [] IIdentifier;
            delete [] AntiIIdentifier;
            
            return SubGrid;
            
        }
        
        
        
        
        void CombineIndices(INT *r_sub,INT *r_int,INT *r,INT *Mask)
        {
            INT ksub=0;
            INT kint=0;
            for(INT i=0; i<this->NCoords; i++)
            {
                if(ksub+kint!=i)
                {
                    std::cerr << "# ERROR IN COMBINING SUB AND INT INDICES" << std::endl;
                    exit(0);
                }
                if(Mask[i]==1)
                {
                    r[i]=r_sub[ksub];
                    ksub++;
                }
                else if(Mask[i]==0)
                {
                    r[i]=r_int[kint];
                    kint++;
                }
            }
        }
        
        
       
        
        
        void IntegralUpdate(INT *ppi_subgrid, complex<DOUBLE> val, complex<DOUBLE> &sum, DOUBLE VolumeElement, INT *IntegralSize, INT IntegralDIM)
        {
            // check if evaluated at boundary //
            DOUBLE weight=1.;
            for(INT i=0; i<IntegralDIM; i++)
            {
                // test if this is a boundary term, if so compute the weight for trapezodial integration //
                if( ppi_subgrid[i]==0 || ppi_subgrid[i]==IntegralSize[i]-1 )
                {
                    weight *= 0.5;
                }
            }
            sum += weight*val*VolumeElement;
        }
        
        
        // ////////////////////////////////////////////////////////////////////////////////////////////////////
        // CREATE A SUBGRID BY TAKING A SLICE AT FIXED COORDINATES
        // /////////////////////////////////////////////////////////////////////////////////////////////////////
        GRID *CreateSubGRID_Slice(INT SliceDIM, INT *Mask, INT *SliceIndx)
        {
            INT NewDim = NCoords-SliceDIM;          // new dimension of subgrid //
            INT *GridSize = new INT[NewDim];
            string *IIdentifier = new string[NewDim];
            
            // initialize GridSize of subgrid here //
            INT k=0;
            for(INT i=0; i<NCoords; i++)
            {
                if(Mask[i]==1)
                {
                    GridSize[k]=this->Size[i];
                    k++;
                }
            }
            
            // create that new subgrid grid //
            GRID *SubGrid = new GRID(NewDim, GridSize, IIdentifier);
            
            // copy coordinates //
            INT subgrid_k=0;
            for(INT i=0; i<NCoords; i++)
            {
                if(Mask[i]==1)
                {
                    SubGrid->Size[subgrid_k] = Size[i];
                    SubGrid->Identifier[subgrid_k] = Identifier[i];
                    
                    for(INT j=0; j<SubGrid->Size[subgrid_k]; j++)
                    {
                        SubGrid->p[subgrid_k][j] = p[i][j];
                    }
                    
                    subgrid_k++;
                }
            }
            
            // copy data //
            INT *ppi = new INT[NCoords];
            INT *ppi_subgrid = new INT[NewDim];
            for(INT i=0; i<this->dataSize; i++)
            {
                this->Coords(i, ppi);
                
                // test if that point is also found on the sub grid //
                INT counter=0;
                for(INT j=0; j<NCoords; j++)
                {
                    if(Mask[j]==1)
                    {
                        // we keep that dimension, so the point is always found on the subgrid //
                        counter++;
                    }
                    else if(Mask[j]==0 && SliceIndx[j]==ppi[j])
                    {
                        // if that dimension is dropped, keep only if at the slice //
                        counter++;
                    }
                }
                
                if(counter==NCoords)
                {
                    // determine point on subgrid from point on grid //
                    INT k=0;
                    for(INT j=0; j<NCoords; j++)
                    {
                            if(Mask[j]==1)
                            {
                                ppi_subgrid[k]=ppi[j];
                                k++;
                            }
                    }
                    
                    // keep this point and sort it into the subgrid //
                    complex<DOUBLE> val = this->Get(ppi);
                    SubGrid->Indx(ppi_subgrid);
                    SubGrid->Set(ppi_subgrid, val);
                }
            }
            
            
            return SubGrid;
        }
        
        
        
        void ShowRanges()
        {
            // Prints the data ranges in each coordinate //
            for(INT i=0; i<NCoords; i++)
            {
                std::cout << "# " <<  Identifier[i] << " ";
                for(INT j=0; j<Size[i]; j++)
                {
                    std::cout << p[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }
        
        
        void Set(INT *pp, complex<DOUBLE> val)
        {
            F[ Indx(pp) ] = val;
        }

        void Set(INT p0, INT p1, INT p2, INT p3, INT p4, INT p5, complex<DOUBLE> val)
        {
            F[ Indx(p0,p1,p2,p3,p4,p5) ] = val;
        }
        
        complex<DOUBLE> Get(INT *pp)
        {
            return F[ Indx(pp) ];
        }

        complex<DOUBLE> Get(INT p0, INT p1, INT p2, INT p3, INT p4, INT p5)
        {
            return F[ Indx(p0,p1,p2,p3,p4,p5) ];
        }
        
        // ////////////////////////////////////////////
        // INDEXING //
        // ////////////////////////////////////////////
        INT Coords(INT Indx, INT *pp)
        {
            INT vol = dataSize;
            INT site = Indx;
            for(INT i=NCoords-1; i>=0; i--)
            {
                vol /= Size[i];
                pp[i] = site / vol;
                site -=  vol*pp[i];
            }
            if(site != 0)
            {
                std::cerr << "de-indexing went wrong";
                exit(0);
            }
        }
        
        INT IntegralCoords(INT IntegralDIM, INT INtegralIndx, INT *Integralpp, INT totalIntegralSize, INT *IntegralSize)
        {
            INT vol = totalIntegralSize;
            INT site = INtegralIndx;
            for(INT i=IntegralDIM-1; i>=0; i--)
            {
                vol /= IntegralSize[i];
                Integralpp[i] = site / vol;
                site -=  vol*Integralpp[i];
            }
            if(site != 0)
            {
                std::cerr << "de-indexing went wrong [INTEGRAL COORDS]";
                exit(0);
            }
        }
        
        
        INT Indx(INT *pp)
        {
            INT result=0;
            INT numerator=1;
            for(INT i=NCoords-1; i>=0; i--)
            {
                numerator *= Size[i];
                result += pp[i]*dataSize/numerator;
            }
            return result;
        }
        
        INT Indx(INT p0, INT p1, INT p2, INT p3, INT p4, INT p5)
        {
            return p5*(dataSize/Size[5]) + p4*(dataSize/(Size[5]*Size[4])) + p3*(dataSize/(Size[5]*Size[4]*Size[3]))
            + p2*(dataSize/(Size[5]*Size[4]*Size[3]*Size[2])) + p1*(dataSize/(Size[5]*Size[4]*Size[3]*Size[2]*Size[1]))
            + p0*(dataSize/(Size[5]*Size[4]*Size[3]*Size[2]*Size[1]*Size[0]));
        }
        // ///////////////////////////////////////////////
        // END OF INDEXING //
        // ////////////////////////////////////////////////
        
        // data //
        INT NCoords;
        INT dataSize;
        INT *Size;
        string *Identifier;
        INT **pInd;
        DOUBLE **p;
        complex<DOUBLE> *F;
        INT *MPI_lower;
        INT *MPI_upper;

        //
    private:
    };
    
}
















//
//
//
//
//
//
//        // INTEGRATES OUT ANY DIMENSION YOU LIKE, SPECIFIED BY BITWISE MASK //
//        // IF Mask[dim]=0, this dimension will be integrated out,
//        // IF Mask[dim]=1, this dimension will be kept
//        // IntegralDIM is the number of dimensions that are to be integrated out //
//        GRID *CreateSubGRID_Integrate(INT IntegralDIM, INT *Mask)
//        {
//            INT NewDim = NCoords-IntegralDIM;          // new dimension of subgrid //
//            INT *GridSize = new INT[NewDim];
//            string *IIdentifier = new string[NewDim];
//
//            // initialize GridSize of subgrid here //
//            INT k=0;
//            for(INT i=0; i<NCoords; i++)
//            {
//                if(Mask[i]==1)
//                {
//                    GridSize[k]=this->Size[i];
//                    k++;
//                }
//            }
//
//            // create that new subgrid grid //
//            GRID *SubGrid = new GRID(NewDim, GridSize, IIdentifier);
//
//            // copy coordinates //
//            INT subgrid_k=0;
//            for(INT i=0; i<NCoords; i++)
//            {
//                if(Mask[i]==1)
//                {
//                    SubGrid->Size[subgrid_k] = Size[i];
//                    SubGrid->Identifier[subgrid_k] = Identifier[i];
//
//                    for(INT j=0; j<SubGrid->Size[subgrid_k]; j++)
//                    {
//                        SubGrid->p[subgrid_k][j] = p[i][j];
//                    }
//                    subgrid_k++;
//                }
//            }
//
//            // Integrate data //
//            INT *ppi = new INT[NCoords];
//            INT *ppi_subgrid = new INT[NewDim];
//
//            INT *IntegralSize = new INT[IntegralDIM];
//            INT totalIntegralSize=1;
//            DOUBLE VolumeElement=1.0;
//            k=0;
//            for(INT i=0; i<NCoords; i++)
//            {
//                if(Mask[i]==0)
//                {
//                    IntegralSize[k] = this->Size[i];
//                    totalIntegralSize *= IntegralSize[k];
//
//                    VolumeElement *= (this->p[i][1]-this->p[i][0]);
//                    k++;
//                }
//            }
//
//            std::cout << "# total integral size " << totalIntegralSize << std::endl;
//            std::cout << "# volume element is " << VolumeElement << std::endl;
//
//            INT *ppi_integral = new INT[IntegralDIM];
//            for(INT i=0; i<SubGrid->dataSize; i++)
//            {
//                SubGrid->Coords(i,ppi_subgrid);
//                DOUBLE sum=0.0;
//
//                for(INT j=0; j<totalIntegralSize; j++)
//                {
//                    // determine the coordinates that are to be integrated over //
//                    IntegralCoords(IntegralDIM, j, ppi_integral, totalIntegralSize, IntegralSize);
//
//                    // store them into some large coordinate //
//                    INT counter_external_coord=0; INT counter_integral_coord=0;
//                    for(INT i=0; i<NCoords; i++)
//                    {
//                        if(Mask[i]==0)
//                        {
//                            ppi[i]=ppi_integral[counter_integral_coord];
//                            counter_integral_coord++;
//                        }
//                        else if(Mask[i]==1)
//                        {
//                            ppi[i]=ppi_subgrid[counter_external_coord];
//                            counter_external_coord++;
//                        }
//                    }
//
//                    // using this complete set of coordinates, read out the value of the integrand //
//                    DOUBLE val = this->Get(ppi);
//
//                    // feed into integration routine //
//                    IntegralUpdate(ppi_subgrid, val, sum, VolumeElement, IntegralSize, IntegralDIM);
//                }
//
//                SubGrid->Set(ppi_subgrid, sum);
//            }
//
//            // CLEAN UP //
////            delete [] GridSize;
////            delete [] IIdentifier;
////            delete [] ppi_integral;
////            delete [] ppi;
////            delete [] ppi_subgrid;
////            delete [] IntegralSize;
//
//            return SubGrid;
//            }
#endif /* data_grids_h */
