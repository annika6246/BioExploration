#include "include/geneML.h"
#include "include/SSTAnalysis.h"
#include "include/sequencePhylogeny.h"


int main() {
    runClassifier("input_files/antibiotic_gene_data.csv");
    SSTAnalysis("input_files/noaa_oni_values.csv");
    processSequences("input_files/myoglobin.fasta");
    return 0;
}