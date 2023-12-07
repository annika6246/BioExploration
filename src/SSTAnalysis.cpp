#include "../include/SSTAnalysis.h"

using namespace mlpack;

void SSTAnalysis(const std::string& filename) {
    // load data
    arma::mat data;
    data::Load(filename, data, true);

    // initialize clusters and run model
    const size_t clusters = 6;
    arma::Row<size_t> assignments;
    KMeans<> kmeans;
    kmeans.Cluster(data, clusters, assignments);

    // output results to a csv file
    std::ofstream output("output_files/kmeans_assignments.csv");
    for (size_t i = 0; i < assignments.n_elem; i++) {
        output << i << "," << assignments[i] << std::endl;
    }
}