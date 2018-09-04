#ifndef random_seed_h
#define random_seed_h

#include<gsl/gsl_rng.h>



//
//namespace SEEDS
//{
//    INT *GenerateLocalSeeds(INT MPI_ID, INT SEEDS_PER_MPI_THREAD, INT LOCAL_MPI_SEED)
//    {
//        INT *MySEED = new INT[SEEDS_PER_MPI_THREAD];
//        
//        LOCAL_MPI_SEED += MPI_ID*429;
//        long int ranLRmax = 99999999999;
//        
//        Random *random;
//        random = new Random();
//        random->gslRandomInit(LOCAL_MPI_SEED);
//        
////        gsl_rng *r;
////        r = gsl_rng_alloc(gsl_rng_mt19937);
////        gsl_rng_set(r, LOCAL_MPI_SEED);
//
//        
//        for(INT i=0; i<SEEDS_PER_MPI_THREAD; i++)
//        {
//            MySEED[i]=gsl_rng_uniform_int(r, ranLRmax);
//        }
//  
//        gsl_rng_free(r);
//    }
//    
//    void PrintSeeds(INT MPI_ID, INT MPI_THREADS, INT SEEDS_PER_MPI_THREAD, INT *SEED_LIST)
//    {
//        for(INT i=0; i<MPI_THREADS; i++)
//        {
//            if(i==MPI_ID)
//            {
//                for(INT j=0; j<SEEDS_PER_MPI_THREAD; j++)
//                {
//                    std::cout << "THREAD " << MPI_ID << " SEED " <<  SEED_LIST[j] << std::endl;
//                }
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//            
//    }
//    
//    
//    
//    
//}






//
//        const gsl_rng_type * rngtype = gsl_rng_default;
//        gsl_rng_env_setup();
//        generator_rng = gsl_rng_alloc(rngtype);
//        gsl_rng_set(global_rng, LOCAL_MPI_SEED);
//





































#endif
