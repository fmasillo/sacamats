
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "sacamats.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
    //std::cerr << argc << '\n';
    if (argc != 2 && argc != 4) {
        //fprintf(stderr, "usage: %s {file of filenames} {isMismatchingSymbolNeeded (1 - yes, 0 - no)} {output file} {verbose}\n", argv[0]);
        fprintf(stderr, "usage: %s {file of filenames}\n or %s {file of filenames} -o {output.txt}\n or ", argv[0], argv[0]);
        exit(1);
    }
    if(argc == 4 && strcmp(argv[2], "-o")){
        std::cout << argv[2] << '\n';
        std::cout << "The only flag that can be added is -o for outputting the GSA.\n";
        exit(1);
    }
    // double totalCollectionSize = 0;
//    double ttotal = getTime();

//    char output_opt = argv[3][0];
    //std::cerr << argv[4] << "\n";
    //bool verbose = argv[4][0] == '1';
    //std::cerr << "verbose = " << verbose << "\n";
    FILE *infilesfile = fopen(argv[1], "r");
    if (!infilesfile) {
        fprintf(stderr, "Error opening file of filenames %s\n", argv[1]);
        exit(1);
    }

    char *filename = (char *) malloc(1024);
    if (!(fgets(filename, 1024, infilesfile))) {
        fprintf(stderr, "Error reading first filename from file of filenames.\n");
        exit(1);
    }

    filename[strlen(filename) - 1] = 0;
    char *refFileName = new char[1024];
    strcpy(refFileName, filename);

    //fprintf(stderr, "About to read SA of ref\n");
    //read the suffix array
    // char safilename[256];
    // sprintf(safilename, "%s.sa", filename);
    // std::cerr << safilename << '\n';
    // infile = fopen(safilename, "r");
    // if (!infile) {
    //     fprintf(stderr, "Error opening suffix array of input file %s\n", safilename);
    //     exit(1);
    // }
    // unsigned int san = 0;
    // fseek(infile, 0, SEEK_END);
    // san = ftell(infile) / sizeof(unsigned int);
    // std::cerr << "san = " << san << '\n';
    // fseek(infile, 0, SEEK_SET);
    // unsigned int *sa = new unsigned int[san];
    // if (san != fread(sa, sizeof(unsigned int), san, infile)) {
    //     fprintf(stderr, "Error reading sa from file\n");
    //     exit(1);
    // }
    // fclose(infile);

    //std::cerr << "Checking SA\n";
    //uint64_t lcpsum = 0;
    //for(uint32_t i = 1; i < san; i++){
    //   if((i%10000000) == 0) std::cerr << "Checked: " << i << '\n';
    //   uint32_t s1 = sa[i-1];
    //   uint32_t s2 = sa[i];
    //   //std::cerr << s1 << " " << s2 << '\n'; 
    //   for(uint32_t j=0;j<n;j++){
    //      //std::cerr << x[s1+j] << " " << x[s2+j] << '\n'; 
    //      if(x[s1+j] != x[s2+j]){
    //         if(x[s1+j] > x[s2+j]){
    //            std::cerr << "Suffixes out of order: " << i << '\n';
    //            exit(1);
    //         }
    //         lcpsum+=j;
    //         break;
    //      }
    //   }
    //}
    //std::cerr << "Avg. LCP: " << lcpsum/n << "\n";

    //fprintf(stderr, "Reference (size = %u):\n\t", n);
    //for (int bi = 0; bi < n; bi++) {
    //    fprintf(stderr, "%u, ", x[bi]);
    //}
    fprintf(stderr, "\n");

    //compute relative LZ factorization
    char * _ = fgets(filename, 1024, infilesfile);
    filename[strlen(filename) - 1] = '\0';
    
    //lzInitialize(x, n, std::stoul(argv[2]), refFileName, filename);
    //lzFactorize(filename, 0, argv[3], verbose);
    std::vector<std::pair<uint32_t, int32_t>> MSGSA;
    computeGSA(refFileName, filename, MSGSA);
    if(argc == 4){
        std::cout << "Saving to file\n";
        //std::cerr << MSGSA.size() << '\n';
        std::ofstream ofp(argv[3], std::ios_base::binary);
        ofp.write(reinterpret_cast<const char*>(&MSGSA[0]), sizeof(std::pair<uint32_t, int32_t>)*MSGSA.size());
        ofp.close();
    }
    //lzInitialize(x, n, refFileName, filename);
    //lzFactorize(filename, 0);

//     int seqno = 1;
//     unsigned int totalNumFactors = 0;
// //    double totalBitsOut = 0;
//     while (fgets(filename, 1024, infilesfile)) {
//         filename[strlen(filename) - 1] = '\0';
//         fprintf(stderr, "Reading sequence %d (%s)\n", seqno, filename);
        
//         fprintf(stderr, "---Factorizing sequence %d\n", seqno);
//         double tlz = getTime();
//         totalNumFactors += lzFactorize(filename, seqno, argv[3], verbose);
//         fprintf(stderr, "---Time to lz factorize sequence %d = %.2f secs\n", seqno, getTime() - tlz);
//         fprintf(stdout, "%.2f,%d,%u,", getTime() - tlz, totalNumFactors, n);
// //        totalBitsOut += bits_out;
//         seqno++;
        
//     }
//    fprintf(stderr, "----Collection parsed into %u factors.\n", totalNumFactors);
//    fprintf(stderr, "----Total time: %.2f secs.\n", getTime() - ttotal);
//    fprintf(stderr, "----Estimated encoded size (minus ref): %.2f megabytes.\n",
//            (totalNumFactors * (ceil(log(n) / log(2))) +
//             totalNumFactors * ceil(log(totalCollectionSize / totalNumFactors) / log(2))) / 8 / (1 << 20));
//    fprintf(stderr, "----Estimated encoded size (incl. ref): %.2f megabytes.\n",
//            (2 * n + totalNumFactors * (ceil(log(n) / log(2))) +
//             totalNumFactors * ceil(log(totalCollectionSize / totalNumFactors) / log(2))) / 8 / (1 << 20));
//    fprintf(stderr, "----Size using Elias Delta (minus ref): %.2f megabytes.\n", (totalBitsOut) / 8 / (1 << 20));
//    fprintf(stderr, "----Size using Elias Delta (incl. ref): %.2f megabytes.\n",
//            (2 * n + totalBitsOut) / 8 / (1 << 20));
    exit(0);

    free(filename);
    //delete[] x;
    //delete[] sa;
    return 0;
}
