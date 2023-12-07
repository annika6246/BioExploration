#include "../include/sequencePhylogeny.h"

using namespace seqan;

void processSequences(const std::string& inputFile) {
    // define sequence type (amino acid for input data)
    typedef String<AminoAcid> TSequence;

    // process data in fasta format
    SeqFileIn seqFileIn(inputFile.c_str());

    // write columns in output csv
    std::ofstream outputCSV("output_files/output_scores.csv");
    outputCSV << "Sequence1,Sequence2,AlignmentScore\n";

    // initialize variables to store species labels and sequence data
    StringSet<CharString> labels;
    StringSet<TSequence> sequences;

    // read data from input file into containers
    readRecords(labels, sequences, seqFileIn);

    // iterate through species
    for (size_t i = 0; i < length(sequences); ++i)
    {
        // for all other species, compute pairwise alignment score, avoids duplicates
        for (size_t j = i + 1; j < length(sequences); ++j)
        {
            // create align instance
            Align<TSequence> align;

            // resize for pairwise alignment
            resize(rows(align), 2);
            assignSource(row(align, 0), sequences[i]);
            assignSource(row(align, 1), sequences[j]);

            // compute alignment score
            int score = globalAlignment(align, Score<int, Simple>(1, -1, -1));

            // output score to csv file
            outputCSV << labels[i] << "," << labels[j] << "," << score << "\n";
        }
    }

    // close files
    close(seqFileIn);
    outputCSV.close();
}