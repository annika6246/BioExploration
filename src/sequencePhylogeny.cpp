#include "../include/sequencePhylogeny.h"

using namespace seqan;

void processSequences(const std::string& inputFile) {
    typedef String<AminoAcid> TSequence;

    SeqFileIn seqFileIn(inputFile.c_str());

    std::ofstream csvFile("output_files/output_scores.csv");
    csvFile << "Sequence1,Sequence2,AlignmentScore\n";

    StringSet<CharString> labels;
    StringSet<TSequence> sequences;

    readRecords(labels, sequences, seqFileIn);

    for (size_t i = 0; i < length(sequences); ++i)
    {
        for (size_t j = i + 1; j < length(sequences); ++j)
        {
            Align<TSequence> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), sequences[i]);
            assignSource(row(align, 1), sequences[j]);

            int score = globalAlignment(align, Score<int, Simple>(1, -1, -1));

            csvFile << labels[i] << "," << labels[j] << "," << score << "\n";
        }
    }

    close(seqFileIn);
    csvFile.close();
}