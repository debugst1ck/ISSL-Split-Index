#include <iostream>

#include <vector>
#include <list>
#include <unordered_map>

#include <algorithm>
#include <execution>
#include <future>
#include <mutex>

#include <string>

#include <fstream>
#include <filesystem>

#include "phmap.h"
#include "cfdPenalties.h"

#include <windows.h>
#include <psapi.h>

// Compile time constants ->

constexpr const char NEWLINE = '\n';

// Custom data structures used ->

/// @brief Enum representing nucleotide types
enum Nucleotide
{
    ADENINE = 0b00,
    CYTOSINE = 0b01,
    GUANINE = 0b10,
    THYMINE = 0b11,
};

enum ScoreMethod
{
    UNKNOWN = 0,
    MIT = 1,
    CFD = 2,
    MIT_AND_CFD = 3,
    MIT_OR_CFD = 4,
    AVG_MIT_CFD = 5,
};

/// @brief Header structure for sub-indices
/// @attention Designed to be fixed 80 bytes long (10 std::uint64_t)
struct SubIndexHeader
{
    std::uint64_t sliceIndex;                             // Current slice number
    std::uint64_t sequenceLength;                         // Length of the sequence
    std::uint64_t sliceWidth;                             // Width of the slice
    std::uint64_t offtargetCount;                         // Number of unique offtargets
    std::uint64_t scoreCount;                             // Total number of precalculated scores
    std::uint64_t sliceCount;                             // Total number of slices
    const std::uint64_t futureMetadata[4] = {0, 0, 0, 0}; // Placeholder metadata, reduce the array number and add member to stuct
};

/// @brief Checks if a file exists at the given path.
/// @param filePath The path of the file to check.
void checkFileExists(const std::filesystem::path &filePath)
{
    if (!std::filesystem::exists(filePath))
    {
        std::cerr << ("Failed to open file: " + filePath.string()) << NEWLINE;
        std::exit(EXIT_FAILURE);
    }
}

/// @attention Changes to the order of the Nucleotide enum values will affect the mapping to characters
/// @brief Converts a Nucleotide enum value to its corresponding character representation
/// @param nucleotide The Nucleotide enum value
/// @return The character representation of the nucleotide
char nucleotideToCharacter(const Nucleotide &nucleotide)
{
    constexpr static const char SIGNATURE_INDEX[4] = {'A', 'C', 'G', 'T'};
    return SIGNATURE_INDEX[nucleotide];
}

/// @brief Converts a character to its corresponding Nucleotide enum value
/// @param character The character representing a nucleotide (A, C, G, T)
/// @return The Nucleotide enum value corresponding to the character
Nucleotide characterToNucleotide(const char &character)
{
    switch (character)
    {
    case 'C':
        return Nucleotide::CYTOSINE;
    case 'G':
        return Nucleotide::GUANINE;
    case 'T':
        return Nucleotide::THYMINE;
    default: // Case 'A' or unwanted stuff should be 0b00
        return Nucleotide::ADENINE;
    }
}

/// @brief Converts a sequence signature to a DNA sequence string
/// @param signature The sequence signature to convert
/// @param sequenceLength The length of the DNA sequence
/// @return DNA sequence string
std::string signatureToSequence(const std::uint64_t &signature, const std::uint64_t &sequenceLength)
{
    std::string sequence;
    sequence.reserve(sequenceLength);
    for (std::uint64_t i = 0; i < sequenceLength; i++)
    {
        sequence += nucleotideToCharacter(static_cast<Nucleotide>((signature >> (i << 1)) & 0b11));
    }
    return sequence;
}

/// @brief Converts the given string of sequence to 2-bit binary encoded signature
/// @param charIterator Iterator to the string, will work in any direction
/// @param sequenceLength Length of the sequence/string
/// @return The sequence encoded in binary
std::uint64_t sequenceToSignature(std::string::const_iterator charIterator, const std::uint64_t &sequenceLength)
{
    std::uint64_t signature = 0;
    for (std::uint64_t i = 0; i < sequenceLength; i++)
    {
        signature |= static_cast<std::uint64_t>(characterToNucleotide(*charIterator++)) << (i << 1);
    }
    return signature;
}

/// @brief This function reads 64-bit integers from a binary file starting from a specified position and stores them in a vector.
/// @param result A reference to a vector of integers to be read from the file will be stored in this vector.
/// @param filePath The path to the binary file to be read.
/// @param position The position in the file to start reading from.
/// @param count Number of 64-bit integers to be read from the file.
/// @return This function returns true if the file was successfully read.
bool read64Bits(std::vector<std::uint64_t> &result, const std::filesystem::path &filePath, const std::uint64_t &position, const std::uint64_t &count)
{
    std::ifstream inputFile(filePath, std::ios::binary);
    inputFile.seekg(position);
    result.resize(count);
    for (int i = 0; i < count; i++)
    {
        inputFile.read(reinterpret_cast<char *>(&result[i]), sizeof(std::uint64_t));
    }
    return true;
}

/// @brief Get specified number of scores from the ISSL sub index file and store them into the parallel hash map
/// @param scoresTable Reference to the parallel hash map for the storage of the mismatches/scores
/// @param fileName Constant reference to path of ISSL file
/// @param scoresCount Number of mismatches/scores to be extracted
/// @param position The offset to the scores table within the ISSL
/// @return `true` if the scores are successfully extracted
bool readScoresTable(phmap::flat_hash_map<std::uint64_t, double> &scoresTable, const std::filesystem::path &fileName, const std::uint64_t &scoresCount, const std::uint64_t &position)
{
    scoresTable.reserve(scoresCount);
    std::ifstream inputFile(fileName, std::ios::binary);
    inputFile.seekg(position);
    char buffer[16];
    for (std::uint64_t scoreIndex = 0; scoreIndex < scoresCount; scoreIndex++)
    {
        inputFile.read(buffer, sizeof(buffer));
        scoresTable.insert(*reinterpret_cast<std::pair<std::uint64_t, double> *>(&buffer));
    }
    return true;
}

/// @brief Get the precalculated header from the ISSL file
/// @param filePath The path to the file
/// @return The header, duh
SubIndexHeader getSubIndexHeader(const std::filesystem::path &filePath)
{
    std::vector<std::uint64_t> headerArray;
    read64Bits(headerArray, filePath, 0, sizeof(SubIndexHeader) >> 3);
    return *reinterpret_cast<SubIndexHeader *>(headerArray.data());
}

/// @brief Parses the string argument provided to enum values for concrete definition
/// @param argument The argument as a string reference
/// @return The `ScoreMethod` enum value.
ScoreMethod parseScoreMethod(const std::string &argument)
{
    static const std::unordered_map<std::string, ScoreMethod> lookUpTable = {
        {"mit", ScoreMethod::MIT},
        {"cfd", ScoreMethod::CFD},
        {"and", ScoreMethod::MIT_AND_CFD},
        {"or", ScoreMethod::MIT_OR_CFD},
        {"avg", ScoreMethod::AVG_MIT_CFD},
    };
    std::unordered_map<std::string, ScoreMethod>::const_iterator searchResult = lookUpTable.find(argument);
    if (searchResult == lookUpTable.end()) // If the argument is not recognized by argument
    {
        return ScoreMethod::UNKNOWN;
    }
    else
    {
        return searchResult->second;
    }
}

/// @brief Determines whether MIT scores should be calculated based on the given ScoreMethod.
/// @param method ScoreMethod Enum
/// @return True if MIT scores should be calculated, false otherwise.
bool getCalculateMIT(const ScoreMethod &method)
{
    switch (method)
    {
    case ScoreMethod::MIT:
        return true;
    case ScoreMethod::MIT_AND_CFD:
        return true;
    case ScoreMethod::MIT_OR_CFD:
        return true;
    case ScoreMethod::AVG_MIT_CFD:
        return true;
    default:
        return false;
    }
}

/// @brief Determines whether CFD scores should be calculated based on the given ScoreMethod.
/// @param method ScoreMethod Enum
/// @return True if CFD scores should be calculated, false otherwise.
bool getCalculateCFD(const ScoreMethod &method)
{
    switch (method)
    {
    case ScoreMethod::CFD:
        return true;
    case ScoreMethod::MIT_AND_CFD:
        return true;
    case ScoreMethod::MIT_OR_CFD:
        return true;
    case ScoreMethod::AVG_MIT_CFD:
        return true;
    default:
        return false;
    }
}

const std::uint64_t computeBitmask(const std::uint64_t &width)
{
    return (1ULL << width) - 1;
}

const std::uint64_t computeLimit(const std::uint64_t &width)
{
    return 1ULL << width;
}

/// @brief Reads DNA sequence signatures from a binary ISSL file and converts them into 64-bit integers.
/// @param querySignatures Reference to a vector to store the read query signatures.
/// @param filePath The path to the binary file containing DNA sequences.
/// @param sequenceLength The length of the DNA sequences.
/// @return The number of query signatures read and stored in the vector.
std::uint64_t readQuerySignatures(
    std::vector<std::uint64_t> &querySignatures,
    const std::filesystem::path &filePath,
    const std::uint64_t &sequenceLength)
{
    std::ifstream inputFile(filePath, std::ios::binary);
    std::string buffer;
    buffer.reserve(sequenceLength);
    while (inputFile >> buffer)
    {
        querySignatures.push_back(sequenceToSignature(buffer.cbegin(), sequenceLength));
    }
    return querySignatures.size();
}

/// @brief Computes the number of mismatches between two DNA signatures.
/// @param searchSignature The signature of the search DNA sequence.
/// @param offtarget The signature of the offtarget DNA sequence.
/// @return The number of mismatches between the search and offtarget DNA sequences.
std::uint64_t computeMismatches(const std::uint64_t searchSignature, const std::uint64_t offtarget)
{
    /** Find the positions of mismatches
     *
     *  Search signature (SS):    A  A  T  T    G  C  A  T
     *                           00 00 11 11   10 01 00 11
     *
     *        Off-target (OT):    A  T  A  T    C  G  A  T
     *                           00 11 00 11   01 10 00 11
     *
     *                SS ^ OT:   00 00 11 11   10 01 00 11
     *                         ^ 00 11 00 11   01 10 00 11
     *                  (XORd) = 00 11 11 00   11 11 00 00
     *
     *        XORd & evenBits:   00 11 11 00   11 11 00 00
     *                         & 10 10 10 10   10 10 10 10
     *                   (eX)  = 00 10 10 00   10 10 00 00
     *
     *         XORd & oddBits:   00 11 11 00   11 11 00 00
     *                         & 01 01 01 01   01 01 01 01
     *                   (oX)  = 00 01 01 00   01 01 00 00
     *
     *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
     *                         | 00 01 01 00   01 01 00 00
     *            mismatches   = 00 01 01 00   01 01 00 00
     *
     *   popcount(mismatches):   4
     */
    const std::uint64_t XORSignatures = searchSignature ^ offtarget;
    const std::uint64_t evenBits = XORSignatures & 0xAAAAAAAAAAAAAAAAULL;
    const std::uint64_t oddBits = XORSignatures & 0x5555555555555555ULL;
    return (evenBits >> 1) | oddBits;
}

std::uint64_t computeSearchSlice(const std::uint64_t search, const std::uint64_t sliceWidth, const std::uint64_t sliceIndex)
{
    const std::uint64_t sliceShift = sliceWidth * sliceIndex;
    const std::uint64_t sliceMask = computeBitmask(sliceWidth);
    return (search >> sliceShift) & sliceMask;
}

///@brief Calculates the CFD scores for a given offtarget based on the distance between the search signature and the offtarget.
///@param distance The number of mismatches between the search signature and offtarget.
///@param maximumDistance The maximum allowed distance for CFD calculation.
///@param searchSignature The signature of the search DNA sequence.
///@param signatureId The identifier for the offtarget signature.
///@param occurrences The number of occurrences of the offtarget.
///@param offtargets Container containing the offtarget signatures.
///@return The CFD score for the given offtarget.
double calculateCFDScores(
    std::uint64_t distance, const std::uint64_t maximumDistance,
    std::uint64_t searchSignature, std::uint32_t signatureId, std::uint32_t occurrences,
    std::vector<std::uint64_t> &offtargets)
{
    double cfdScore = 0.0;
    // If there are no mismatches, set CFD score to 1
    if (distance == 0)
    {
        cfdScore = 1;
    }
    // If there are mismatches within the allowed range, calculate CFD score
    else if (distance > 0 && distance <= maximumDistance)
    {
        // Initialize CFD score with PAM penalty for NGG (TODO: avoid hard-coding PAM)
        cfdScore = cfdPamPenalties[0b1010];
        // Iterate over each position in the DNA sequence
        for (size_t position = 0; position < 20; position++)
        {
            // Extract nucleotide identity at the current position from the search signature
            std::uint64_t searchSigIdentityPosition = (searchSignature & (0b11UL << (position * 2))) >> (position * 2);
            searchSigIdentityPosition <<= 2;
            // Extract nucleotide identity at the current position from the offtarget
            std::uint64_t offtargetIdentityPos = (offtargets[signatureId] & (0b11UL << (position * 2))) >> (position * 2);
            // Check if the nucleotide identity differs between search and offtarget
            if (searchSigIdentityPosition >> 2 != offtargetIdentityPos)
            {
                // Generate a mask combining position, search signature identity, and offtarget identity
                std::uint64_t mask = (position << 4 | searchSigIdentityPosition | (offtargetIdentityPos ^ 0b11UL));
                // Multiply the CFD score by the penalty associated with the mask
                cfdScore *= cfdPosPenalties[mask];
            }
        }
    }
    return cfdScore * (double)occurrences;
}

/// @brief Prints a progress bar to the console to represent the progress of a task.
/// @param progress The current progress value.
/// @param total The total value representing the completion of the task.
/// @param barWidth The width of the progress bar (default is 50).
void printProgressBar(int progress, int total, int barWidth = 50)
{
    float percentage = static_cast<float>(progress) / static_cast<float>(total);
    int progressBarWidth = static_cast<int>(percentage * barWidth);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < progressBarWidth)
            std::cout << "=";
        else
            std::cout << " ";
    }

    std::cout << "] " << std::setw(5) << (percentage * 100) << "%\r";
    std::cout.flush();
}

/// @brief Reads the offsets to slices from a ISSL file and stores them in a vector.
/// @param offsetToSlice Reference to a vector to store the slice offsets.
/// @param filePath The path to the binary file.
/// @param initialPosition The initial position in the file to start reading from.
/// @param sliceLimit The limit on the number of slices to read.
/// @return The position in the file after reading the slice offsets.
std::streampos readSliceOffsets(
    std::vector<std::uint64_t> &offsetToSlice,
    const std::filesystem::path &filePath,
    const std::streampos &initialPosition,
    const std::uint64_t &sliceLimit)
{
    std::ifstream inputFile(filePath, std::ios::binary);
    inputFile.seekg(initialPosition);
    for (std::uint64_t slice = 0x00; slice < sliceLimit; slice++)
    {
        std::uint64_t sliceCount = 0;
        offsetToSlice.push_back(inputFile.tellg());
        inputFile.read(reinterpret_cast<char *>(&sliceCount), sizeof(std::uint64_t));
        inputFile.seekg(static_cast<std::uint64_t>(inputFile.tellg()) + (sizeof(std::uint64_t) * sliceCount));
    }
    return inputFile.tellg();
}

/// @brief Scores a slice based on given parameters and updates scores for MIT and CFD.
/// @param isslFile The path to the binary ISSL file.
/// @param offset The offset in the ISSL file to read the slice from.
/// @param searchSignature The signature of the search DNA sequence.
/// @param maximumDistance The maximum allowed distance for scoring.
/// @param calculateCFD Flag indicating whether to calculate CFD scores.
/// @param calculateMIT Flag indicating whether to calculate MIT scores.
/// @param offtargets Reference to a vector containing offtarget DNA signatures.
/// @param scoresTable Hash map containing precomputed scores for MIT.
/// @param offtargetToggles Vector indicating whether an offtarget has been processed to avoid duplicate scoring.
/// @return A pair of scores (MIT score, CFD score) for the given slice.
std::pair<double, double> scoreSlice(
    const std::filesystem::path isslFile, const std::streampos offset,
    const std::uint64_t searchSignature, const std::uint64_t maximumDistance,
    const bool calculateCFD, const bool calculateMIT,
    std::vector<std::uint64_t> &offtargets,
    phmap::flat_hash_map<std::uint64_t, double> &scoresTable,
    std::vector<bool> &offtargetToggles)
{
    std::ifstream inputFile(isslFile, std::ios::binary);
    inputFile.seekg(offset);
    std::uint64_t sliceKeyCount = 0;
    inputFile.read(reinterpret_cast<char *>(&sliceKeyCount), sizeof(std::uint64_t));
    double scoreMit = 0.0;
    double scoreCfd = 0.0;
    for (std::uint64_t keyIndex = 0; keyIndex < sliceKeyCount; keyIndex++)
    {
        std::uint64_t key;
        inputFile.read(reinterpret_cast<char *>(&key), sizeof(std::uint64_t));
        std::uint32_t signatureId = key & 0xFFFFFFFFULL;
        std::uint32_t occurrences = (key >> (32));
        std::uint64_t mismatches = computeMismatches(searchSignature, offtargets[signatureId]);
        std::uint64_t distance = std::__popcount(mismatches);

        if (distance >= 0 && distance <= maximumDistance)
        {
            if (!offtargetToggles[signatureId])
            {
                if (calculateMIT && (distance > 0))
                {
                    scoreMit += scoresTable[mismatches] * static_cast<double>(occurrences);
                }
                if (calculateCFD)
                {
                    scoreCfd += calculateCFDScores(distance, maximumDistance, searchSignature, signatureId, occurrences, offtargets);
                }
                offtargetToggles[signatureId] = true;
            }
        }
    }
    return std::make_pair(scoreMit, scoreCfd);
}

/// @brief Asynchronously processes a search query using sub-indices and computes MIT and CFD scores for the respective query sig.
/// @param searchSignature The query signature of the search DNA sequence.
/// @param subIndexHeaders Vector of sub-index headers containing slice information.
/// @param subIndexPaths Vector of file paths to sub-indices.
/// @param offsetsToSlices Vector of vectors containing offsets to slices in each sub-index.
/// @param maximumDistance The maximum allowed distance for scoring.
/// @param calculateCFD Flag indicating whether to calculate CFD scores.
/// @param calculateMIT Flag indicating whether to calculate MIT scores.
/// @param offtargets Reference to a vector containing offtarget DNA signatures.
/// @param scoresTable Hash map containing precomputed scores for MIT.
/// @return A pair of scores (MIT score, CFD score) for the given query.
std::pair<double, double> processQuery(
    const std::uint64_t &searchSignature,
    const std::vector<SubIndexHeader> &subIndexHeaders,
    const std::vector<std::filesystem::path> &subIndexPaths,
    const std::vector<std::vector<std::uint64_t>> &offsetsToSlices,
    const std::uint64_t maximumDistance,
    const bool calculateCFD, const bool calculateMIT,
    std::vector<std::uint64_t> &offtargets,
    phmap::flat_hash_map<std::uint64_t, double> &scoresTable)
{
    std::vector<std::future<std::pair<double, double>>> sliceScores(subIndexHeaders[0].sliceCount);
    std::vector<bool> offtargetToggles(subIndexHeaders[0].offtargetCount, false);
    for (std::uint64_t sliceIndex = 0; sliceIndex < subIndexHeaders[0].sliceCount; sliceIndex++)
    {
        std::uint64_t searchSlice = computeSearchSlice(searchSignature, subIndexHeaders[sliceIndex].sliceWidth, sliceIndex);
        sliceScores[sliceIndex] = std::async(std::launch::async, scoreSlice, subIndexPaths[sliceIndex], offsetsToSlices[sliceIndex][searchSlice], searchSignature, maximumDistance, calculateCFD, calculateMIT, std::ref(offtargets), std::ref(scoresTable), std::ref(offtargetToggles));
    }
    double totalScoreMit = 0.0;
    double totalScoreCfd = 0.0;
    for (std::uint64_t sliceIndex = 0; sliceIndex < subIndexHeaders[0].sliceCount; sliceIndex++)
    {
        std::pair<double, double> score = sliceScores[sliceIndex].get();
        totalScoreMit += score.first;
        totalScoreCfd += score.second;
    }
    return std::make_pair((10000.0 / (100.0 + totalScoreMit)), (10000.0 / (100.0 + totalScoreCfd)));
}

int main(int argc, char **argv)
{
    const bool debug = true;
    const std::chrono::time_point start = std::chrono::steady_clock::now(); // Steady clock is a monotonic clock

    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " [CANDIDATE_GUIDE_PATH] [MAXIMUM_DISTANCE] [SCORE_THRESHOLD] [ISSL_TABLE]" << NEWLINE;
        return EXIT_FAILURE;
    }

    const std::filesystem::path queryFile = std::filesystem::path(argv[1]);
    const std::uint64_t maximumDistance = std::stoi(argv[2]);
    const double threshold = std::stod(argv[3]);
    const std::string scoreMethodArg = std::string(argv[4]);
    std::vector<std::filesystem::path> subIndexPaths;

    for (int i = 5; i < argc; i++)
    {
        subIndexPaths.push_back(std::filesystem::path(argv[i]));
    }

    const ScoreMethod scoreMethod = parseScoreMethod(scoreMethodArg);

    const bool calculateCFD = getCalculateCFD(scoreMethod);
    const bool calculateMIT = getCalculateMIT(scoreMethod);

    std::vector<SubIndexHeader> subIndexHeaders;

    for (std::filesystem::path subIndexPath : subIndexPaths)
    {
        subIndexHeaders.push_back(getSubIndexHeader(subIndexPath));
    }

    // Choose a file to read in offtargets and scores
    std::vector<std::uint64_t> offtargets;
    read64Bits(offtargets, subIndexPaths[0], sizeof(SubIndexHeader), subIndexHeaders[0].offtargetCount);

    const std::uint64_t offsetToScores = sizeof(SubIndexHeader) + subIndexHeaders[0].offtargetCount * 8;
    phmap::flat_hash_map<std::uint64_t, double> scoresTable;
    readScoresTable(scoresTable, subIndexPaths[0], subIndexHeaders[0].scoreCount, offsetToScores);

    std::vector<std::uint64_t> querySignatures;
    readQuerySignatures(querySignatures, queryFile, subIndexHeaders[0].sequenceLength);

    const std::uint64_t offsetToSliceList = (sizeof(SubIndexHeader) + (subIndexHeaders[0].offtargetCount * 8) + (subIndexHeaders[0].scoreCount * 16));

    std::vector<std::vector<std::uint64_t>> offsetsToSlices(subIndexHeaders[0].sliceCount);

    for (std::uint64_t i = 0; i < subIndexHeaders[0].sliceCount; i++)
    {
        readSliceOffsets(offsetsToSlices[i], subIndexPaths[i], offsetToSliceList, computeLimit(subIndexHeaders[0].sliceWidth));
    }

    std::vector<std::pair<double, double>> scores(querySignatures.size());

    for (std::uint64_t searchIndex = 0; searchIndex < querySignatures.size(); searchIndex++)
    {
        scores[searchIndex] = processQuery(std::ref(querySignatures[searchIndex]), std::ref(subIndexHeaders), std::ref(subIndexPaths), std::ref(offsetsToSlices), maximumDistance, calculateCFD, calculateMIT, std::ref(offtargets), std::ref(scoresTable));
        printProgressBar(searchIndex + 1, querySignatures.size(), 100);
    }
    std::cout << NEWLINE; // End progress bar

    constexpr const char fieldSeparator[3] = ", ";
    for (std::uint64_t searchIndex = 0; searchIndex < querySignatures.size(); searchIndex++)
    {
        std::string querySequence = signatureToSequence(querySignatures[searchIndex], subIndexHeaders[0].sequenceLength);
        std::pair<double, double> score = scores[searchIndex];
        double scoreMIT = calculateMIT ? score.first : -1.0;
        double scoreCFD = calculateCFD ? score.second : -1.0;
        std::cout << querySequence << fieldSeparator << scoreMIT << fieldSeparator << scoreCFD << NEWLINE;
    }

    // If "--debug" parameter set, memory tracking only for windows sadly
    if (debug)
    {
        const std::chrono::microseconds runtime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        PROCESS_MEMORY_COUNTERS_EX pmc;
        GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc));

        std::cerr << querySignatures.size() << fieldSeparator;
        std::cerr << subIndexHeaders[0].offtargetCount << fieldSeparator;
        std::cerr << subIndexHeaders[0].sliceWidth << fieldSeparator;
        std::cerr << (static_cast<double>(pmc.PeakWorkingSetSize) / (1024 * 1024)) << " MB" << fieldSeparator; // Converting B to MB
        std::cerr << runtime.count() << " microseconds" << NEWLINE;
    }
    return EXIT_SUCCESS;
}
