
#ifndef BIOEXPLORATION_GENEML_H
#define BIOEXPLORATION_GENEML_H

#include <core/data/load.hpp>
#include <methods/random_forest/random_forest.hpp>
#include <core/data/confusion_matrix.hpp>
#include <core/data/confusion_matrix_impl.hpp>
#include <core/cv/k_fold_cv.hpp>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

void readCSV(const std::string& filename, std::vector<int>& labels, std::vector<std::string>& geneData);
std::vector<std::vector<int>> processSequences(std::vector<std::string>& geneData);

void matrixToImage(const arma::Mat<size_t>& cM);
void runClassifier(const std::string& filename);

#endif //BIOEXPLORATION_GENEML_H
