/// @file CreateSplitIndex.cpp
/// @attention Requires the POSIX threads library for Linux environments, add the -pthreads argument to g++ if not compiling

#include <iostream>

#include <string>

#include <list>
#include <vector>
#include <queue>
#include <map>

#include <cstdint>
#include <limits>

#include <filesystem>
#include <fstream>

#include <thread>

#include <cstdlib> // For error codes

#include <chrono> // For measuring execution times

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

// Error catching functions ->

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

// Utility Functions ->

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

/// @attention Changes to the order of the Nucleotide enum values will affect the mapping to characters
/// @brief Converts a Nucleotide enum value to its corresponding character representation
/// @param nucleotide The Nucleotide enum value
/// @return The character representation of the nucleotide
char nucleotideToCharacter(const Nucleotide &nucleotide)
{
    constexpr static const char SIGNATURE_INDEX[4] = {'A', 'C', 'G', 'T'};
    return SIGNATURE_INDEX[nucleotide];
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

const std::uint64_t computeBitmask(const std::uint64_t &width)
{
    return (1ULL << width) - 1;
}

const std::uint64_t computeLimit(const std::uint64_t &width)
{
    return 1ULL << width;
}

/// @brief The number of slice lists required is calculated
/// @param sequenceLength The length of sequence/offtarget
/// @param width The sequence width in bits (2,4,8,16...)
/// @return The number of slices
const std::uint64_t computeSliceCount(const std::uint64_t &sequenceLength, const std::uint64_t &width)
{
    return (sequenceLength << 1) / width;
}

// Scoring functions, I have no idea how to document these ones

/// @brief Calculates the local score for a sequence with mismatch patterns
/// @param mismatchArray A vector containing the positions of mismatches in the sequence
/// @param length The length of the sequence
/// @return The calculated local score
double calculateLocalScore(const std::vector<std::uint64_t> &mismatchArray, const std::uint64_t &length)
{
    double T1 = 1.0, T2, T3, d = 0.0, score;
    /* Mismatch penalty array */
    double M[] = {0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

    /* 1st term */
    for (std::uint64_t i = 0; i < length; ++i)
    {
        T1 = T1 * (1.0 - M[mismatchArray[i]]);
    }

    /* 2nd term */
    if (length == 1)
        d = 19.0;
    else
    {
        for (std::uint64_t i = 0; i < length - 1; ++i)
        {
            d += mismatchArray[i + 1] - mismatchArray[i];
        }
        d = d / (length - 1);
    }
    T2 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1);

    /* 3rd term */
    T3 = 1.0 / (length * length);

    /* Local score */
    score = T1 * T2 * T3 * 100;
    return score;
}

/// @brief Predicts the local score for a sequence with a given XOR of signatures
/// @param XORSignatures The XOR of signatures representing the mismatch pattern
/// @param sequenceLength The length of the DNA sequence
/// @return The predicted local score for the sequence
double predictLocalScore(std::uint64_t XORSignatures, std::uint64_t &sequenceLength)
{
    std::vector<std::uint64_t> mismatchArray(20);
    std::uint64_t m = 0;
    for (std::uint64_t j = 0; j < sequenceLength; j++)
    {
        if ((XORSignatures >> (j << 1)) & 0b11)
        {
            mismatchArray[m++] = j;
        }
    }
    if (m == 0)
    {
        return 0.0;
    }
    return calculateLocalScore(mismatchArray, m);
}

/// @brief Computes the mismatch masks for a sequence of a given length with a specified number of mismatches
/// @param seqLength The length of the sequence
/// @param mismatches The number of mismatches to generate masks for
/// @return A vector containing mismatch masks encoded in 2-bit format
std::vector<uint64_t> computeMasksTwoBit(int seqLength, int mismatches)
{
    std::vector<uint64_t> masks;
    if (mismatches < seqLength)
    {
        if (mismatches > 0)
        {
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches - 1))
            {
                masks.push_back((1LLU << (seqLength - 1) * 2) + mask);
            }
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches))
            {
                masks.push_back(mask);
            }
        }
        else
        {
            masks.push_back(0LLU);
        }
    }
    else
    {
        uint64_t tempMask = 0;
        for (int i = 0; i < seqLength; i++)
        {
            tempMask |= (1LLU << i * 2);
        }
        masks.push_back(tempMask);
    }
    return masks;
}

// File interactions

/// @brief Get the size of a line in a file, including line termination characters
/// @param filePath The path to the input file
/// @param sequenceLength The length of the sequence in the line
/// @return The size of the line, including line termination characters
std::uint64_t computeLineSize(const std::filesystem::path &filePath, const std::uint64_t &sequenceLength)
{
    std::ifstream inputFile(filePath, std::ios::binary);
    std::string line;
    std::getline(inputFile, line);
    return line.back() == '\r' ? sequenceLength + 2 : sequenceLength + 1; // +2 for \r\n and +1 for \n
}

/// @brief Parses the given arguments to filepaths
/// @attention Does not check for file extension as ".issl" as of now
/// @param start Start index of the arguments
/// @param end Number of arguments
/// @param pathList The raw arguments for the program
/// @return The list of sub-index paths as doubly-linked list
std::list<std::filesystem::path> generateSubIndexPaths(int start, int end, char **pathList)
{
    std::list<std::filesystem::path> subIndexFilePaths;
    for (std::uint64_t i = start; i < end; i++)
    {
        std::filesystem::path file = std::filesystem::path(pathList[i]);
        std::ofstream temp(file, std::ios::trunc);
        subIndexFilePaths.emplace_back(file);
    }
    return subIndexFilePaths;
}

// Generative functions

/// @brief Reads in the "offtarget" file sequentially, convert them to binary and counts its occurrences
/// @attention This function assumes that offtargets file is sorted, otherwise invalid data
/// @param offtargetsTable A empty table of offtargets encoded in 2-bit binary, and the occurrences of the offtargets, to be filled
/// @param filePath The file-path to the offtargets file in raw unicode format
/// @param sequenceLength The length of the sequence/offtarget
/// @param sequenceCount The estimated number of sequences in the file (don't really have to be accurate)
/// @return The number of unique offtargets in the file
const std::uint64_t constructOfftargetsTable(
    std::vector<std::pair<std::uint64_t, std::uint64_t>> &offtargetsTable,
    const std::filesystem::path &filePath,
    const std::uint64_t &sequenceLength, const std::uint64_t &sequenceCount)
{
    offtargetsTable.reserve(sequenceCount);                                    // Preallocate memory from sequenceCount for worst case scenario.
    std::uint64_t uniqueSiteIndex = std::numeric_limits<std::uint64_t>::max(); // Vector index, initialized max (0-1) for overflow to 0 for first index
    std::ifstream inputFile(filePath, std::ios::binary);
    std::uint64_t signature, previousSignature;
    std::string buffer;
    while (inputFile >> buffer)
    {
        signature = sequenceToSignature(buffer.cbegin(), sequenceLength);
        if (previousSignature == signature)
        {
            offtargetsTable[uniqueSiteIndex].second++;
        }
        else
        {
            uniqueSiteIndex++;
            offtargetsTable[uniqueSiteIndex] = {signature, 1};
        }
        previousSignature = signature;
    }
    return uniqueSiteIndex + 1; // 0 index to count
}

/// @brief Constructs the sub-indices by assigning signatures ID and occurrences to their respective slice lists
/// @attention The structure of sliceList is [sliceIndex][sliceValue][sequenceSignatureKey]
/// @param sliceLists A vector of vector of vector of std::uint64_t representing the slice lists
/// @param sequenceTableIterator An iterator pointing to the sequence table
/// @param offtargetsCount The number of unique offtargets
/// @param sliceWidth The width of each slice in bits
/// @param sliceCount The total number of slice lists (equal to the number of sub-index files)
void constructIndex(
    std::vector<std::vector<std::vector<std::uint64_t>>> &sliceLists,
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator sequenceTableIterator,
    const std::uint64_t &offtargetsCount, const std::uint64_t &sliceWidth, const std::uint64_t &sliceCount)
{
    const std::uint64_t bitmask = computeBitmask(sliceWidth);
    for (std::uint64_t signatureIndex = 0; signatureIndex < offtargetsCount; signatureIndex++)
    {
        const std::uint64_t signature = sequenceTableIterator->first;
        const std::uint64_t occurrence = sequenceTableIterator->second;
        for (std::uint64_t sliceIndex = 0; sliceIndex < sliceCount; sliceIndex++)
        {
            const std::uint64_t sliceValue = (signature >> (sliceIndex * sliceWidth)) & bitmask;
            const std::uint64_t sequenceSignatureKey = (occurrence << 32) | signatureIndex;
            sliceLists[sliceIndex][sliceValue].push_back(sequenceSignatureKey);
        }
        sequenceTableIterator++;
    }
}

/// @brief Constructs a table of pre-calculated local scores for different mismatch patterns
/// @param scoresTable A map containing pre-calculated local scores indexed by their corresponding mismatch pattern
/// @param sequenceLength The length of the DNA sequence for which scores are being calculated
/// @param sliceWidth The width of each slice in bits, used to generate mismatch patterns
/// @return The total number of pre-calculated scores
std::uint64_t constructScoresTable(
    std::map<std::uint64_t, double> &scoresTable,
    std::uint64_t sequenceLength, std::uint64_t sliceWidth)
{
    std::uint64_t maxDist = sequenceLength * 2 / sliceWidth - 1;
    std::uint64_t scoresCount = 0;
    for (int i = 1; i <= maxDist; i++)
    {
        std::vector<std::uint64_t> tempMasks;
        tempMasks = computeMasksTwoBit(20, i);
        for (auto mask : tempMasks)
        {
            double score = predictLocalScore(mask, sequenceLength);
            scoresTable[mask] = score;
            scoresCount++;
        }
    }
    return scoresCount;
}

// Main writing functions

/// @brief Writes the offtarget sequence signatures to the specified file
/// @param filePath The path to the output file where offtarget sequence signatures will be written
/// @param sequenceTableIterator An iterator pointing to the sequence table, (signatures and occurrences)
/// @param offtargetCount The total number of unique offtarget sequence signatures
void writeOfftargets(
    const std::filesystem::path &filePath,
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator sequenceTableIterator,
    const std::uint64_t &offtargetCount)
{
    std::ofstream outputFile(filePath, std::ios::binary | std::ios::app);
    for (std::uint64_t i = 0; i < offtargetCount; i++)
    {
        outputFile.write(reinterpret_cast<const char *>(&(sequenceTableIterator->first)), sizeof(std::uint64_t));
        sequenceTableIterator++;
    }
}

/// @brief Writes the sequence signature keys for a specific slice value to the specified file
/// @attention Writes the size of the slice first then the contents (ID and occurrences) of the slice
/// @param filePath The path to the output file where the sequence signature keys will be written
/// @param sliceListIterator An iterator pointing to the slice list containing sequence signature keys
/// @param sliceLimit The maximum number of slices in a single slice list. 2^x
void writeSliceList(
    const std::filesystem::path &filePath,
    std::vector<std::vector<std::uint64_t>>::const_iterator sliceListIterator,
    const std::uint64_t &sliceLimit)
{
    std::ofstream outputFile(filePath, std::ios::binary | std::ios::app);
    for (std::uint64_t signatureIndex = 0; signatureIndex < sliceLimit; signatureIndex++)
    {
        std::uint64_t sliceSize = sliceListIterator->size();
        outputFile.write(reinterpret_cast<const char *>(&sliceSize), sizeof(std::uint64_t));
        for (std::uint64_t sequenceSignatureKey : *sliceListIterator)
        {
            outputFile.write(reinterpret_cast<const char *>(&sequenceSignatureKey), sizeof(std::uint64_t));
        }
        sliceListIterator++;
    }
}

/// @brief Writes the pre-calculated scores to the specified file
/// @param filePath The path to the output file where pre-calculated scores will be written
/// @param scoresTable A map containing pre-calculated scores indexed by their corresponding mismatch pattern
void writeScores(
    const std::filesystem::path &filePath,
    std::map<std::uint64_t, double> &scoresTable)
{
    std::ofstream outputFile(filePath, std::ios::binary | std::ios::app);
    for (std::pair<const std::uint64_t, double> mask_score : scoresTable)
    {
        outputFile.write(reinterpret_cast<const char *>(&mask_score.first), sizeof(std::uint64_t));
        outputFile.write(reinterpret_cast<const char *>(&mask_score.second), sizeof(double));
    }
}

/// @brief Writes sub-index headers to the specified sub-index files
/// @param filePathsIterator An iterator pointing to the list of sub-index file paths
/// @param headerIterator An iterator pointing to the list of sub-index headers
/// @param sliceCount The total number of sub-index files
void writeHeaders(
    std::list<std::filesystem::path>::const_iterator filePathsIterator,
    std::list<SubIndexHeader>::const_iterator headerIterator,
    const std::uint64_t &sliceCount)
{
    for (std::uint64_t i = 0; i < sliceCount; i++)
    {
        std::fstream outputFile(*filePathsIterator, std::ios::binary | std::ios::in | std::ios::out);
        outputFile.write(reinterpret_cast<const char *>(&*headerIterator), sizeof(SubIndexHeader));
        headerIterator++;
        filePathsIterator++;
    }
}

// Here the multithreading functions implemented

/// @brief Stops the main thread temporarily for all the worker threads to catch up
/// @param workers The queue of threads to synchronize with main thread
void synchronize(std::queue<std::thread> &workers)
{
    while (!workers.empty()) // TODO: Prevent infinite loop if possible
    {
        if (workers.front().joinable())
        {
            workers.front().join();
            workers.pop();
        }
    }
}

/// @brief Utilizes the std::thread to create threads and assigns them the task to write the offtargets
/// @details Used iterators for using the list in any order.
/// @param filePathsIterator Iterator to the list of files
/// @param offtargetsTable Iterator to table of offtargets and occurrences (only reads the offtargets signature)
/// @param sliceCount Number of sub-index files
/// @param offtargetCount Number of offtargets
/// @param workers The thread worker is moved into the workers queue for access
void concurrentlyWriteOfftargets(
    std::list<std::filesystem::path>::const_iterator filePathsIterator,
    const std::vector<std::pair<std::uint64_t, std::uint64_t>>::const_iterator offtargetsTableIterator,
    const std::uint64_t &sliceCount, const std::uint64_t &offtargetCount,
    std::queue<std::thread> &workers)
{
    for (std::uint64_t i = 0; i < sliceCount; i++)
    {
        std::thread worker(writeOfftargets, *filePathsIterator, offtargetsTableIterator, offtargetCount);
        workers.emplace(std::move(worker));
        filePathsIterator++;
    }
}

/// @brief Concurrently writes the pre-calculated scores to multiple sub-index files using multiple threads
/// @param filePathsIterator An iterator pointing to the list of sub-index file paths
/// @param scoresTable A map containing pre-calculated scores indexed by their corresponding mismatch pattern
/// @param sliceCount The total number of sub-index files
/// @param workers A queue of thread objects where worker threads will be added for synchronization
void concurrentlyWriteScores(
    std::list<std::filesystem::path>::const_iterator filePathsIterator,
    std::map<std::uint64_t, double> &scoresTable,
    const std::uint64_t &sliceCount,
    std::queue<std::thread> &workers)
{
    for (std::uint64_t i = 0; i < sliceCount; i++)
    {
        std::thread worker(writeScores, *filePathsIterator, std::ref(scoresTable));
        workers.emplace(std::move(worker));
        filePathsIterator++;
    }
}

/// @brief Concurrently writes the slice lists to multiple sub-index files using multiple threads
/// @param filePathsIterator An iterator pointing to the list of sub-index file paths
/// @param sliceLists A 3D vector containing the slice lists of sequence signature keys
/// @param sliceCount The total number of sub-index files
/// @param sliceWidth The width of each slice in bits
/// @param workers A queue of thread objects where worker threads will be added for synchronization
void concurrentlyWriteSliceLists(
    std::list<std::filesystem::path>::const_iterator filePathsIterator,
    const std::vector<std::vector<std::vector<std::uint64_t>>> &sliceLists,
    const std::uint64_t &sliceCount, const std::uint64_t &sliceWidth,
    std::queue<std::thread> &workers)
{
    for (std::uint64_t i = 0; i < sliceCount; i++)
    {
        std::thread worker(writeSliceList, *filePathsIterator, sliceLists[i].cbegin(), computeLimit(sliceWidth));
        workers.emplace(std::move(worker));
        filePathsIterator++;
    }
}

/// @brief The main program. Duh
/// @param argc Number of arguments
/// @param argv The arguments as pointer to strings
/// @return exit code
int main(int argc, char **argv)
{
    const std::chrono::time_point start = std::chrono::steady_clock::now(); // Steady clock is a monotonic clock
    // Start content here

    // If arguments insufficient
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0] << " [offtargetSites.txt] [sequence length] [slice width (bits)] [slice paths]" << NEWLINE;
        return EXIT_FAILURE;
    }

    // Reading data in
    const std::filesystem::path filePath = std::filesystem::path(argv[1]);
    const std::uint64_t sequenceLength = std::stoi(argv[2]);
    const std::uint64_t sliceWidth = std::stoi(argv[3]);

    // Preliminary check
    checkFileExists(filePath);

    const std::uint64_t lineSize = computeLineSize(filePath, sequenceLength);
    const std::uint64_t sliceCount = computeSliceCount(sequenceLength, sliceWidth);
    const std::uint64_t fileSize = std::filesystem::file_size(filePath);
    const std::uint64_t sequenceCount = fileSize / lineSize;
    const std::uint64_t sliceLimit = computeLimit(sliceWidth);

    // Secondary checks
    if (sliceCount != (argc - 4))
    {
        std::cerr << "Number of output files insufficient. Required: " << sliceCount << ", Given: " << (argc - 4) << NEWLINE;
        return EXIT_FAILURE;
    }
    if (fileSize % lineSize != 0)
    {
        std::cerr << "File does is not a multiple of the expected line length " << lineSize << NEWLINE;
        return EXIT_FAILURE;
    }
    if (sequenceLength > 32)
    {
        std::cerr << "Sequence length is greater than 32, which is the maximum supported currently" << NEWLINE;
        return EXIT_FAILURE;
    }

    // Create dem files from arguments
    std::cout << "Wiping and creating given files" << NEWLINE;
    std::list<std::filesystem::path> subIndexPaths = generateSubIndexPaths(4, argc, argv);

    std::cout << "Creating headers" << NEWLINE;
    std::list<SubIndexHeader> headerList;
    for (std::uint64_t i = 0; i < sliceCount; i++)
    {
        SubIndexHeader header;
        header.sliceIndex = i;
        header.sequenceLength = sequenceLength;
        header.sliceWidth = sliceWidth;
        header.sliceCount = sliceCount;
        headerList.emplace_back(header);
    }

    // thread pool
    std::queue<std::thread> workers;

    std::cout << "Writing empty headers" << NEWLINE;
    writeHeaders(subIndexPaths.cbegin(), headerList.cbegin(), sliceCount);

    std::cout << "Reading, generating and writing offtargets" << NEWLINE;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> offtargetsTable;
    const std::uint64_t offtargetsCount = constructOfftargetsTable(offtargetsTable, filePath, sequenceLength, sequenceCount);
    concurrentlyWriteOfftargets(subIndexPaths.cbegin(), offtargetsTable.cbegin(), sliceCount, offtargetsCount, workers);

    std::cout << "Generating and writing pre-calculated scores" << NEWLINE;
    std::map<std::uint64_t, double> scoresTable;
    const std::uint64_t scoresCount = constructScoresTable(scoresTable, sequenceLength, sliceWidth);
    synchronize(workers);
    concurrentlyWriteScores(subIndexPaths.cbegin(), scoresTable, sliceCount, workers);

    std::cout << "Generating and writing slice lists to slice files" << NEWLINE;
    std::vector<std::vector<std::vector<std::uint64_t>>> sliceLists(sliceCount, std::vector<std::vector<std::uint64_t>>(sliceLimit, std::vector<std::uint64_t>()));
    constructIndex(sliceLists, offtargetsTable.cbegin(), offtargetsCount, sliceWidth, sliceCount);
    synchronize(workers);
    concurrentlyWriteSliceLists(subIndexPaths.cbegin(), sliceLists, sliceCount, sliceWidth, workers);

    std::cout << "Updating header with new data and editing it into file" << NEWLINE;
    for (std::list<SubIndexHeader>::iterator header = headerList.begin(); header != headerList.end(); header++)
    {
        header->offtargetCount = offtargetsCount;
        header->scoreCount = scoresCount;
    }
    synchronize(workers);
    writeHeaders(subIndexPaths.cbegin(), headerList.cbegin(), sliceCount);

    // Benchmarking part
    const std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Time taken by program: " << duration.count() << " microseconds" << NEWLINE;
    return EXIT_SUCCESS;
}