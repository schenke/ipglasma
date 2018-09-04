

#ifndef linear_interoplation_h
#define linear_interoplation_h




namespace LINEAR_INTERPOLATION
{
    class Interpolator
    {
    public:
        Interpolator(INT NNCORDS)
        {
            NCORDS=NNCORDS;
            HYPERCUBE=new INT[NCORDS];
            HYPERCUBE_TEMP=new INT[NCORDS];
            
            NTERMS=1;
            for(INT i=0;i<NCORDS;i++){NTERMS*=2;}
            
            BITMASK=new INT*[NTERMS];
            for(INT i=0;i<NTERMS;i++)
            {
                BITMASK[i]=new INT[NCORDS];
                for(INT j=0; j<NCORDS; j++)
                {
                    BITMASK[i][j]=(i>>j)&1;
                }
            }
            
            A = new COMPLEX[NTERMS];
        }
        
        virtual ~Interpolator()
        {
            delete HYPERCUBE;
            delete HYPERCUBE_TEMP;
            delete BITMASK;
        }
        
        
        COMPLEX Get(DOUBLE *x, DATA::GRID *Data)
        {
            FindHypercube(x,Data);
            GetHypercubeValues(A,Data);
            
            
            COMPLEX sum=(0.0,0.0);
            for(INT i=0; i<NTERMS; i++)
            {
                COMPLEX temp=A[i];
                for(INT j=0; j<NCORDS; j++)
                {
                    // determine x2 and x1 //
                    DOUBLE xlower=Data->p[j][HYPERCUBE[j]];
                    DOUBLE xupper=Data->p[j][HYPERCUBE[j]+1];
                    
                    if( BITMASK[i][j]==0 )
                    {
                        temp *= (xupper-x[j])/(xupper-xlower);
                    }
                    else if(BITMASK[i][j]==1)
                    {
                        temp *= (x[j]-xlower)/(xupper-xlower);
                    }
                }
                sum+=temp;
            }
            return sum;
        }
        
        void GetHypercubeValues(COMPLEX *Val, DATA::GRID *Data)
        {
            for(INT i=0;i<NTERMS;i++)
            {
                // COMPUTE INDEX //
                for(INT j=0;j<NCORDS;j++)
                {
                    HYPERCUBE_TEMP[j]=HYPERCUBE[j]+BITMASK[i][j];
                }
                
                // HYPERCUBE_TEMP TO BE USED TO ADDRESS DATA //
                Val[i]=Data->Get(HYPERCUBE_TEMP);
            }
        }
        
        
        
        void FindHypercube(DOUBLE *x, DATA::GRID *Data)
        {
            for(INT i=0;i<NCORDS;i++)
            {
                if( x[i]<Data->p[i][0] || x[i]>Data->p[i][Data->Size[i]-1])
                {
                    std::cerr << i << ": " << x[i] << " - " << Data->p[i][0] << " / " << Data->p[i][Data->Size[i]-1] << std::endl;
                std:cerr << "cannot interpolate, data out of bounds" << std::endl;
                    exit(0);
                }
                
                HYPERCUBE[i]=0;
                INT j=0;
                while( Data->p[i][j] < x[i] )
                {
                    HYPERCUBE[i]=j;
                    j++;
                }
            }
        }
        
        void TestHyperCubeSearch(DOUBLE *x, DATA::GRID *Data)
        {
            FindHypercube(x,Data);
            std::cout << "# data was" << std::endl;
            for(INT i=0; i<NCORDS; i++)
            {
                std::cout << Data->p[i][HYPERCUBE[i]] << " < " <<  x[i] << " < " << Data->p[i][HYPERCUBE[i]+1] << std::endl ;
            }
            
        }
        
        INT NCORDS;
        INT NTERMS;
        INT *HYPERCUBE;
        INT *HYPERCUBE_TEMP;
        INT **BITMASK;
        COMPLEX *A;
    };
    
}











































#endif
