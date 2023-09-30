#include <iostream>

#include <vector>
#include <list>
#include <unordered_map>

#include <string>

#include <fstream>
#include <filesystem>

#include "phmap.h"
#include "cfdPenalties.h"

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

std::uint64_t computeMismatches(const std::uint64_t searchSignature, const std::uint64_t offtarget)
{
    const std::uint64_t XORSignatures = searchSignature ^ offtarget;
    const std::uint64_t evenBits = XORSignatures & 0xAAAAAAAAAAAAAAAAULL;
    const std::uint64_t oddBits = XORSignatures & 0x5555555555555555ull;
    return (evenBits >> 1) | oddBits;
}

std::uint64_t computeSearchSlice(const std::uint64_t search, const std::uint64_t sliceWidth, const std::uint64_t sliceIndex)
{
    const std::uint64_t sliceShift = sliceWidth * sliceIndex;
    const std::uint64_t sliceMask = computeBitmask(sliceWidth);
    return (search >> sliceShift) & sliceMask;
}

double calculateCFDScores(
    std::uint64_t distance, const std::uint64_t maximumDistance,
    std::uint64_t searchSignature, std::uint32_t signatureId, std::uint32_t occurrences,
    std::vector<std::uint64_t> &offtargets)
{
    double cfdScore = 0.0;
    if (distance == 0)
    {
        cfdScore = 1;
    }

    else if (distance > 0 && distance <= maximumDistance)
    {
        cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM
        for (size_t pos = 0; pos < 20; pos++)
        {

            std::uint64_t searchSigIdentityPos = (searchSignature & (0b11UL << (pos * 2))) >> (pos * 2);
            searchSigIdentityPos <<= 2;
            std::uint64_t offtargetIdentityPos = (offtargets[signatureId] & (0b11UL << (pos * 2))) >> (pos * 2);
            if (searchSigIdentityPos >> 2 != offtargetIdentityPos)
            {
                std::uint64_t mask = (pos << 4 | searchSigIdentityPos | (offtargetIdentityPos ^ 0b11UL));
                cfdScore *= cfdPosPenalties[mask];
            }
        }
    }
    return cfdScore * (double)occurrences;
}

bool shouldCheckNextSlice(ScoreMethod scoreMethod, double totalScoreMIT, double totalScoreCFD, double maximumSum)
{
    switch (scoreMethod)
    {
    case ScoreMethod::MIT_AND_CFD:
        return (totalScoreMIT <= maximumSum) && (totalScoreCFD <= maximumSum);

    case ScoreMethod::MIT_OR_CFD:
        return (totalScoreMIT <= maximumSum) || (totalScoreCFD <= maximumSum);

    case ScoreMethod::AVG_MIT_CFD:
        return ((totalScoreMIT + totalScoreCFD) / 2.0) <= maximumSum;

    case ScoreMethod::MIT:
        return totalScoreMIT <= maximumSum;

    case ScoreMethod::CFD:
        return totalScoreCFD <= maximumSum;

    default:
        // Handle an invalid or unknown ScoreMethod here.
        return true;
    }
}

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

std::streampos readSlice(std::vector<std::uint64_t> &slice, const std::filesystem::path &filePath, const std::streampos &initialPosition)
{
    std::ifstream inputFile(filePath, std::ios::binary);
    inputFile.seekg(initialPosition);
    std::uint64_t sliceKeyCount = 0;
    inputFile.read(reinterpret_cast<char *>(&sliceKeyCount), sizeof(std::uint64_t));
    slice.resize(sliceKeyCount);
    for (std::uint64_t keyIndex = 0; keyIndex < sliceKeyCount; keyIndex++)
    {
        inputFile.read(reinterpret_cast<char *>(&slice[keyIndex]), sizeof(std::uint64_t));
    }
    return inputFile.tellg();
}

int main(int argc, char **argv)
{
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

    std::uint64_t offtargetToggleCount = (subIndexHeaders[0].offtargetCount / (sizeof(std::uint64_t) * 8)) + 1;

    std::vector<double> querySignatureMitScores(querySignatures.size());
    std::vector<double> querySignatureCfdScores(querySignatures.size());

    for (std::uint64_t searchIndex = 0; searchIndex < querySignatures.size(); searchIndex++)
    {
        printProgressBar(searchIndex, querySignatures.size(), 100);

        std::vector<std::uint64_t> offtargetToggles(offtargetToggleCount);

        std::uint64_t searchSignature = querySignatures[searchIndex];

        double totalScoreMit = 0.0;
        double totalScoreCfd = 0.0;

        std::uint64_t offTargetSitesScoredCount = 0;

        const double maximumSum = (10000.0 - threshold * 100) / threshold;

        bool checkNextSlice = true;

        for (std::uint64_t sliceIndex = 0; sliceIndex < subIndexHeaders[0].sliceCount; sliceIndex++)
        {
            const SubIndexHeader header = getSubIndexHeader(subIndexPaths[sliceIndex]);
            std::vector<std::uint64_t> slice;

            std::uint64_t searchSlice = computeSearchSlice(searchSignature, header.sliceWidth, sliceIndex);

            readSlice(slice, subIndexPaths[sliceIndex], offsetsToSlices[sliceIndex][searchSlice]);

            for (std::uint64_t keyIndex = 0; keyIndex < slice.size(); keyIndex++)
            {
                std::uint64_t sequenceSignatureKey = slice[keyIndex];

                std::uint32_t signatureId = sequenceSignatureKey & 0xFFFFFFFFULL;
                std::uint32_t occurrences = (sequenceSignatureKey >> (32));

                std::uint64_t mismatches = computeMismatches(searchSignature, offtargets[signatureId]);

                std::uint64_t distance = std::__popcount(mismatches);

                if (distance >= 0 && distance <= maximumDistance)
                {
                    // note that this is a reverse iterator, so addition is reversed
                    std::reverse_iterator<std::vector<std::uint64_t>::iterator> offtargetFlagIterator = offtargetToggles.rbegin() + (signatureId / 64);
                    std::uint64_t seenOfftargetAlready = (*offtargetFlagIterator >> (signatureId % 64)) & 1ULL;

                    if (!seenOfftargetAlready)
                    {
                        if (calculateMIT && (distance > 0))
                        {
                            totalScoreMit += scoresTable[mismatches] * static_cast<double>(occurrences);
                        }
                        if (calculateCFD)
                        {
                            totalScoreCfd += calculateCFDScores(distance, maximumDistance, searchSignature, signatureId, occurrences, offtargets);
                        }

                        *offtargetFlagIterator |= (1ULL << (signatureId % 64));
                        offTargetSitesScoredCount += occurrences;

                        if (!shouldCheckNextSlice(scoreMethod, totalScoreMit, totalScoreCfd, maximumSum))
                        {
                            checkNextSlice = false;
                            break;
                        }
                    }
                }
            }

            if (!checkNextSlice)
            {
                break;
            }
        }
        // TODO: Directly write output to CSV file
        querySignatureMitScores[searchIndex] = 10000.0 / (100.0 + totalScoreMit);
        querySignatureCfdScores[searchIndex] = 10000.0 / (100.0 + totalScoreCfd);
    }

    // Break the progress bar
    std::cout << NEWLINE;
    
    constexpr const char fieldSeparator[] = ", ";
    for (std::uint64_t searchIndex = 0; searchIndex < querySignatures.size(); searchIndex++)
    {
        std::string querySequence = signatureToSequence(querySignatures[searchIndex], subIndexHeaders[0].sequenceLength);
        double scoreMIT = calculateMIT ? querySignatureMitScores[searchIndex] : -1.0;
        double scoreCFD = calculateCFD ? querySignatureCfdScores[searchIndex] : -1.0;
        std::cout << querySequence << fieldSeparator << scoreMIT << fieldSeparator << scoreCFD << NEWLINE;
    }

    const std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cerr << "Time taken by program: " << duration.count() << " microseconds" << NEWLINE;
    return EXIT_SUCCESS;
}
