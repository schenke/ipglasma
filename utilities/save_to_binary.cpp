#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argv, char* argc[])
{
    int Nevent = 50000;
    char inputfile[128];
    sprintf(inputfile, "coords_12C_NCSM_mixed.dat");
    FILE* infile;
    infile = fopen(inputfile,"r");
    char stemp1[100];
    char** stemp2;
    double x, y, z;
    int event_loop_flag = 1;
    int count_event_number = 0;
    // open file for output
    std::string binary_output_filename = "coords_12C_NCSM_mixed.bin";
    remove(binary_output_filename.c_str());
    FILE *outbin = NULL;
    outbin = fopen(binary_output_filename.c_str(), "wb");
    
    for (int iev = 0; iev < Nevent; iev ++) {
        if(feof(infile)) {
            cout << " End the event loop ~~~ " << endl;
            break;
        }
        for (int ii=0; ii < 11; ii++) {
            fscanf(infile,"%lf %lf %lf ", &x, &y, &z);
            float array[] = {
                static_cast<float>(x),
                static_cast<float>(y), 
                static_cast<float>(z),
            };
            fwrite(array, sizeof(float), 3, outbin);
        }
        fscanf(infile,"%lf %lf %lf\n", &x, &y, &z);
        float array[] = {
            static_cast<float>(x),
            static_cast<float>(y), 
            static_cast<float>(z),
        };
        fwrite(array, sizeof(float), 3, outbin);
    }

    fclose(infile);
    fclose(outbin);
    return 0;
}

