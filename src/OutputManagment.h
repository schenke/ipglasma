//
//  OutputManagment.h
//  SubNucleon diffraction
//
//  Created by Niklas Mueller on 3/3/18.
//  Copyright © 2018 Heikki. All rights reserved.
//

#ifndef OutputManagment_h
#define OutputManagment_h

namespace  IO {
    
    std::string OutputDirectory="./";
    std::string OutputFile="test.txt";
    
    template<typename GenericArgument>
    void SetOutputDirectory(GenericArgument x){
        
        OutputDirectory=StringManipulation::StringCast(x,"/");
        
       // if(MPIBasic::ID==0){
        std::cerr << "#OUTPUT DIRECTORY IS " << OutputDirectory <<  std::endl;
       // }
    }
    
    template<typename GenericArgument>
    void SetOutputFile(GenericArgument x){
        
        OutputFile=StringManipulation::StringCast(x);
        
        // if(MPIBasic::ID==0){
        //std::cerr << "#OUTPUT FILE IS " << OutputFile <<  std::endl;
        // }
    }
    
}



#endif /* OutputManagment_h */
