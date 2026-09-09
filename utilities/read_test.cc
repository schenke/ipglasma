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
    
    // Open the binary file
    std::ifstream InStream;
    InStream.precision(15);
    InStream.open(binary_output_filename, std::ios::in | std::ios::binary);

    // Check if the file was successfully opened
    if (!InStream.is_open()) {
        std::cerr << "Error opening file: " << binary_output_filename << std::endl;
        return 2; // Return an error code
    }
    
    for (int iev = 0; iev < Nevent; iev ++) {
        if (InStream.eof())
            break;
        cout << iev << " ";
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 3; j++) {
                float temp;
                InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
                cout << temp << " ";
            }
        }
        cout << endl;
        if (InStream.eof())
            break;
    }

    fclose(infile);
    InStream.close();
    return 0;
}

