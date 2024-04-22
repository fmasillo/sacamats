
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
#include <getopt.h>

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " [options] <input filename>" << std::endl;
  std::cout << "<input filename> is the name of the file containing paths to the reference sequence (in the first line) and to the collection file (in the second line)." << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-p \tread only a prefix of the file expressed in number of characters, def. whole file" << std::endl
        << "\t-t \tnumber of threads to use, def. max available" << std::endl
        << "\t-o \tbasename for the output files, def. <input filename>" << std::endl
        << "\t-h \tprints this help" << std::endl;
  exit(-1);
}

// function for parsing the input arguments
void parseArgs( int argc, char** argv, Args& arg ) {
    int c;
    extern int optind;

    puts("==== Command line:");
    for(int i=0;i<argc;i++)
        printf(" %s",argv[i]);
    puts("");

    std::string sarg;
    while ((c = getopt( argc, argv, "p:t:o:h") ) != -1) {
        switch(c) {
            case 'p':
                sarg.assign(optarg);
                arg.prefixLength = atoll(sarg.c_str()); break;
                // store the prefix length
            case 't':
                sarg.assign(optarg);
                arg.nThreads = atoll(sarg.c_str()); break;
                // store the number of threads
            case 'o':
                sarg.assign(optarg);
                arg.outname.assign(sarg); break;
                // store the output files path
            case 'h':
                print_help(argv); exit(-1);
                // fall through
            default:
                std::cout << "Unknown option. Use -h for help." << std::endl;
                exit(-1); 
        }
    }

    // the only input parameter is the file name
    if (argc == optind+1) {
        arg.filename.assign( argv[optind] );
    }
    else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_help(argv);
    }
    // set output files basename
    if(arg.outname == "") arg.outname = arg.filename + ".gsa";

    std::cout << "==== Parameters:" << std::endl;
    std::cout << "Input file: " << arg.filename << std::endl;
    std::cout << "Output basename: " << arg.outname << std::endl;
    std::cout << "Prefix length: " << arg.prefixLength << std::endl;
    std::cout << "Number of threads: " << arg.nThreads << std::endl;
    std::cout << std::endl;

}

int main(int argc, char **argv) {
    
    // parse the input arguments
    Args arg;
    parseArgs(argc, argv, arg);

    FILE *infilesfile = fopen(arg.filename.c_str(), "r");
    if (!infilesfile) {
        fprintf(stderr, "Error opening file of filenames %s\n", arg.filename.c_str());
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
    fprintf(stderr, "\n");

    //compute relative LZ factorization
    char * _ = fgets(filename, 1024, infilesfile);
    filename[strlen(filename) - 1] = '\0';
    
    std::vector<std::pair<uint32_t, int32_t>> MSGSA;
    computeGSA(refFileName, filename, arg.prefixLength, arg.nThreads, MSGSA);
    
    std::cout << "Saving to file\n";
    //std::cerr << MSGSA.size() << '\n';
    std::ofstream ofp(arg.outname.c_str(), std::ios_base::binary);
    ofp.write(reinterpret_cast<const char*>(&MSGSA[0]), sizeof(std::pair<uint32_t, int32_t>)*MSGSA.size());
    ofp.close();

    free(filename);
    //delete[] x;
    //delete[] sa;
    return 0;
}
