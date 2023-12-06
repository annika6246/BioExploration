cmake_minimum_required(VERSION 3.20)
project(BioExploration)

include_directories("/mnt/c/Program Files (x86)")
include_directories("/mnt/c/Program Files (x86)/mlpack")
include_directories("/mnt/c/Program Files (x86)/mlpack/lib")
include_directories("/mnt/c/Program Files (x86)/mlpack/methods")
include_directories("/mnt/c/Program Files (x86)/mlpack/core")
include_directories("/mnt/c/Program Files (x86)/mlpack/armadillo-12.6.6/armadillo-12.6.6/include")
include_directories("/mnt/c/Program Files (x86)/mlpack/ensmallen-2.21.0/ensmallen-2.21.0/include")
include_directories("/mnt/c/Program Files (x86)/mlpack/cereal-1.3.2/cereal-1.3.2/include")
include_directories("/mnt/c/Program Files (x86)/seqan-library-2.4.0/seqan-library-2.4.0/include")


set(CMAKE_CXX_STANDARD 17)

add_executable(BioExploration driver.cpp
        include/geneML.h
        src/geneML.cpp
        include/SSTAnalysis.h
        src/SSTAnalysis.cpp
        include/sequencePhylogeny.h
        src/sequencePhylogeny.cpp
)