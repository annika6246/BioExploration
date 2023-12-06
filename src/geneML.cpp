#include "../include/geneML.h"

using namespace mlpack;
using namespace arma;

void readCSV(const std::string& filename, std::vector<int>& labels, std::vector<std::string>& geneData) {
    std::ifstream input(filename);
    std::string line{};
    if (!input) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    while (input) {
        std::getline(input, line);
        int ind {};
        ind = line.find("\r", 0);
        if (ind != std::string::npos) {
            line.erase(ind, 2);
        }

        std::string index, label, sequence, temp;

        std::stringstream lineStr(line);

        std::getline(lineStr, index, ',');
        std::getline(lineStr, label, ',');
        std::getline(lineStr, sequence);

        labels.push_back(std::stoi(label));
        geneData.push_back(sequence);
    }
}

std::vector<std::vector<int>> processSequences(std::vector<std::string>& geneData) {
    std::unordered_map<char, int> charMap;
    std::vector<std::vector<int>> encodedGenes;

    double i = 1;
    for (std::string& sequence : geneData) {
        std::vector<int> encodedSequence;
        for (char c : sequence) {
            if (charMap.find(c) == charMap.end()) {
                charMap[c] = i++;
            }
            encodedSequence.push_back(charMap[c]);
        }
        encodedGenes.push_back(encodedSequence);
    }
    return encodedGenes;
}

void runClassifier(const std::string& filename) {
    std::vector<int> oldLabels;
    std::vector<std::string> geneData;

    // read and process data
    readCSV(filename, oldLabels, geneData);
    std::vector<std::vector<int>> processedData = processSequences(geneData);

    const size_t rows = processedData.size();
    const size_t cols = processedData[0].size();

    // cast labels to floats and transpose data
    arma::mat dataset(cols, rows);
    for (size_t i = 0; i < processedData.size(); ++i) {
        for (size_t j = 0; j < processedData[0].size(); ++j) {
            dataset(j, i) = static_cast<double>(processedData[i][j]);
        }
    }

    arma::Row<size_t> labels = arma::conv_to<arma::Row<size_t>>::from(oldLabels);

    // split training and testing set
    arma::mat trainX, testX;
    arma::Row<size_t> trainY, testY;
    data::Split(dataset, labels, trainX, testX, trainY, testY, 0.2);

    // initialize parameters for model
    const size_t numClasses = 2;
    const size_t minimumLeafSize = 2;
    size_t numTrees = 20;
    const double minimumGainSplit = 1e-7;
    size_t maxDepth = 10;

    // initialize random forest model
    RandomForest<GiniGain, RandomDimensionSelect> rf(trainX, trainY, numClasses, numTrees,
                                                     minimumLeafSize, minimumGainSplit, maxDepth);

    // optional hyperparameter tuning
//    HyperParameterTuner<RandomForest<GiniGain, RandomDimensionSelect>, Accuracy, KFoldCV> hpt(5, dataset, labels, numClasses);
//    std::vector<size_t> maxDepthValues = {5, 10, 15};
//    std::vector<size_t> numTreesValues = {10, 20, 50};
//    std::tie(numTrees, maxDepth) = hpt.Optimize(numTreesValues, maxDepthValues);

    const size_t k = 10;
    KFoldCV<RandomForest<GiniGain, RandomDimensionSelect>, Accuracy> cv(k, dataset, labels, numClasses);
    double cvAccuracy = cv.Evaluate(numTrees, minimumLeafSize, minimumGainSplit, maxDepth);

    // evaluate model
    arma::Row<size_t> predictions;
    rf.Classify(testX, predictions);
    double tPrecision = Precision<Binary>::Evaluate(rf, testX, testY);
    double tRecall = Recall<Binary>::Evaluate(rf, testX, testY);
    double tF1 = F1<Binary>::Evaluate(rf, testX, testY);

    // write evaluations to file
    std::ofstream evaluations("Source/bio_files/evaluations.txt");
    if (evaluations.is_open()) {
        evaluations << "10 K-fold cross validated accuracy: " << cvAccuracy << std::endl;
        evaluations << "Precision for testing set: " << tPrecision << std::endl;
        evaluations << "Recall for testing set: " << tRecall << std::endl;
        evaluations << "F1 score for testing set: " << tF1 << std::endl;
        evaluations.close();
    }

    Mat<size_t> output;
    data::ConfusionMatrix(testY, predictions, output, numClasses);
    output.save("output_files/confusion_matrix.csv", csv_ascii);
}