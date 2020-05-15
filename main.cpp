
#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include <unordered_map>

////////////////////////////////////////////////////////////////////////////////////////////////////////Helper functions and objects.
namespace {

//@brief Type used to represent a codon
using Codon = std::array<char, 3>;

char STOP = '.'; //< Stop codons are decoded as the value and not translated
char START = 'M'; //< Start codon is decoded as the value and it is included in the output string
Codon UNKNOWN = {'?', '?', '?'}; //< The value is used as an unknown codon;

//@brief Helper structure that calculates a hash of the Codon type
struct CodonHasher
{
  // Minimum size of size_t requred by c++ standard is 65535, 
  // which is 16 bits of precision unsigned intiger.
  std::size_t operator()(const Codon& cod) const
  {
    //The requested letters can be encoded on 5 bits and I am shifting them to prevent hash colisions
    return std::hash<int>()((cod[0]- 'A')*1024 ^ (cod[1]-'A')*32 + (cod[2]- 'A'));
  }
};

//@brief Helper hash map describing translation from a codon to a protein
std::unordered_map< Codon, char, CodonHasher> codon2prot {
{{'U','U','U'}, {'F'}},   {{'U','U','C'}, {'F'}},   {{'U','U','A'}, {'L'}},   {{'U','U','G'}, {'L'}},   {{'U','C','U'}, {'S'}},   
{{'U','C','C'}, {'S'}},   {{'U','C','A'}, {'S'}},   {{'U','C','G'}, {'S'}},   {{'U','A','U'}, {'Y'}},   {{'U','A','C'}, {'Y'}},   
{{'U','A','A'}, {STOP}},   {{'U','A','G'}, {STOP}},   {{'U','G','U'}, {'C'}},   {{'U','G','C'}, {'C'}},   {{'U','G','A'}, {STOP}},   
{{'U','G','G'}, {'W'}},   {{'C','U','U'}, {'L'}},   {{'C','U','C'}, {'L'}},   {{'C','U','A'}, {'L'}},   {{'C','U','G'}, {'L'}},   
{{'C','C','U'}, {'P'}},   {{'C','C','C'}, {'P'}},   {{'C','C','A'}, {'P'}},   {{'C','C','G'}, {'P'}},   {{'C','A','U'}, {'H'}},  
{{'C','A','C'}, {'H'}},   {{'C','A','A'}, {'Q'}},   {{'C','A','G'}, {'Q'}},   {{'C','G','U'}, {'R'}},   {{'C','G','C'}, {'R'}},   
{{'C','G','A'}, {'R'}},   {{'C','G','G'}, {'R'}},   {{'A','U','U'}, {'I'}},   {{'A','U','C'}, {'I'}},   {{'A','U','A'}, {'I'}},   
{{'A','U','G'}, {'M'}},   {{'A','C','U'}, {'T'}},   {{'A','C','C'}, {'T'}},   {{'A','C','A'}, {'T'}},   {{'A','C','G'}, {'T'}},   
{{'A','A','U'}, {'N'}},   {{'A','A','C'}, {'N'}},   {{'A','A','A'}, {'K'}},   {{'A','A','G'}, {'K'}},   {{'A','G','U'}, {'S'}},   
{{'A','G','C'}, {'S'}},   {{'A','G','A'}, {'R'}},   {{'A','G','G'}, {'R'}},   {{'G','U','U'}, {'V'}},   {{'G','U','C'}, {'V'}},   
{{'G','U','A'}, {'V'}},   {{'G','U','G'}, {'V'}},   {{'G','C','U'}, {'A'}},   {{'G','C','C'}, {'A'}},   {{'G','C','A'}, {'A'}},  
{{'G','C','G'}, {'A'}},   {{'G','A','U'}, {'D'}},   {{'G','A','C'}, {'D'}},   {{'G','A','A'}, {'E'}},   {{'G','A','G'}, {'E'}},   
{{'G','G','U'}, {'G'}},   {{'G','G','C'}, {'G'}},   {{'G','G','A'}, {'G'}},   {{'G','G','G'}, {'G'}}
};

//@brief Helper function to validate one input character
//@returns true if the character is A, U, C or G; false otherwise
bool correct_char(const char& c)
{
  return c == 'A' || c == 'U' || c == 'C' || c == 'G';
}

//@brief Helper function to validate a codon
//@returns true if all its characters are correct
bool is_valid(const Codon& codon)
{
  return correct_char(codon[0]) && correct_char(codon[1]) && correct_char(codon[2]);
}

//@brief Helper function that prints a help message
void print_help()
{
  std::cout << "rosalind_prot" << " [options]" << std::endl <<
      "Options:" << std::endl <<
      "-h | --help           Print this help" << std::endl <<
      "-i | --input          Input file (default in.txt)" << std::endl <<
      "-o | --output         Output file (default out.txt)" << std::endl <<
      "-f | --frame          RNA reading frame. See: https://en.wikipedia.org/wiki/Genetic_code" << std::endl << 
      "                      (chapter: Features, subchapter: Reading frame) Valid: 0, 1, 2 (default: 0)" << std::endl <<
      "-w | --waitforstart   Flag if proteins decoded before the first start codon should be ignored" << std::endl <<
      "                      Valid: 0, 1 (default: 1)" << std:: endl <<
      "-s | --usestop        Flag if proteins decoded after stop codon and before the next start codon" << std::endl << 
      "                      should be ignored. Valid: 0, 1 (default: 1)" << std::endl << 
      std::endl <<
      "The application reads RNA sequence from the input file and writes the decoded proteins to the output file" << std::endl <<
      "Task description: http://rosalind.info/problems/prot/" << std::endl <<
      "The program returns 0 on success an an error code otherwise." << std::endl <<
      std:: endl <<
      "Error code | Meaning " << std::endl <<
      "-1         | Input file path does not have a filename or the file does not exist." << std::endl <<
      "-2         | Output file path does not have a filename." << std::endl <<
      "-3         | Invalid reading frame value" << std::endl <<
      "-4         | Error while openening input and/or output file." << std::endl <<
      "-5         | Input file contains an invalid codon - a codon that contains characters other than A,U,C,G" << std::endl <<
      "-6         | Missing end codon - the decoded sequence in the output file is not finished" << std::endl <<
      "-7         | Error on parsing input parameters" << std::endl;

}

} // namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////

namespace fs = std::experimental::filesystem;

//@brief The most important function here. It decodes the input and writed the decoded values to output
//@param input Valid, open input stream
//@param output Valid, open output stream.
//@param waitForStart Flag if proteins decoded before the first start codon should be ignored
//@param useEndCodon Flag if proteins decoded after stop codon and before the next start codon shpuld be ignored
//@return 0 on success, -5 on invalid codon and -6 if the sentence is not closed by an end codon
int decode(std::ifstream& input, std::ofstream& output, bool waitForStart, bool useEndCodon)
{
  Codon codon = UNKNOWN;
  bool started = !waitForStart;
  bool ended = false;

  // The codon table is reused in every loop iteration. The approach has pros and cons.
  // pros: it does not need to be created and destroyed all over again which is more time-efficient.
  // cons: the physical memory can "wear out" if the loop runs for a long long time, because 
  //       flash memories has limited number of IO operations. It may lead to unexpected reading errors, 
  //       but the issue is unlikely in this application.
  while (input.read(&codon[0], 3)) 
  {
    if(!is_valid(codon))
    {
      if(input.eof())
      {
        continue; // if it's the end of the string
      }
      // return error if it occurs in the middle
      std::cerr << "Invalid codon: " << codon[0] << codon[1] << codon[2] << std::endl;
      return -5;
    }

    char prot = codon2prot.at(codon);

    if(!started)
    {
      if(prot != START)
      {
        continue;
      }
      started = true; // start codon codes also a protein
      ended = false;
    }

    if(prot == STOP)
    {
      if(useEndCodon)
      {
        ended = true;
        started = false;
      }
      continue; // stop codon codes stop only
    }

    if(started && !ended)
    {
      output << prot;
    }

    codon = UNKNOWN;
  }
  
  if(useEndCodon && started && !ended)
  {
    std::cerr << "Missing end codon - the decoded sequence is not finished" << std::endl;
    return -6;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  ////////// default values
  fs::path inFile ("in.txt");
  fs::path outFile("out.txt");
  int offset = 0;
  bool waitForStart = true;
  bool useEndCodon = true;

  ////////// parsing input parameters
  bool customInput = false;
  bool customOutput = false;
  bool customFrame = false;
  bool customStartFlag = false;
  bool customStopFlag = false;
  bool correctParams = true;
  int argNum = 1;
  while(argNum < argc)
  {
    std::string arg(argv[argNum]);
    if((arg == "-i" || arg == "--input") && !customInput && argNum + 1 < argc)
    {
      customInput = true;
      ++argNum;
      inFile = argv[argNum];
      std::cout << "Custom input " << argv[argNum] << std::endl;
    }
    else if((arg == "-o" || arg == "--output") && !customOutput && argNum + 1 < argc)
    {
      customOutput = true;
      ++argNum;
      outFile = argv[argNum];
      std::cout << "Custom output " << argv[argNum] << std::endl;
    }
    else if((arg == "-f" || arg == "--frame") && !customFrame && argNum + 1 < argc)
    {
      customFrame = true;
      ++argNum;
      std::string val = argv[argNum];
      if(val != "0" && val != "1" && val != "2")
      {
        correctParams = false;
        break;
      }
      offset = std::stoi(val);
    }
    else if((arg == "-w" || arg == "--waitforstart") && !customStartFlag && argNum + 1 < argc)
    {
      customStartFlag = true;
      ++argNum;
      std::string val = argv[argNum];
      if(val != "0" && val != "1")
      {
        correctParams = false;
        break;
      }
      waitForStart = (val == "1");
    }
    else if((arg == "-s" || arg == "--usestop") && !customStopFlag && argNum + 1 < argc)
    {
      customStopFlag = true;
      ++argNum;
      std::string val = argv[argNum];
      if(val != "0" && val != "1")
      {
        correctParams = false;
        break;
      }
      useEndCodon = (val == "1");
    }
    else if((arg == "-h" || arg == "--help"))
    {
      print_help();
    }
    else
    {
      correctParams = false;
      break;
    }
    ++argNum;
  }

  if(!correctParams)
  {
    std::cout << "Error parsing input arguments: " << std::endl;
    print_help();
    return -7;
  }

  ////////// parameters validation 
  if(!inFile.has_filename() || !fs::exists(inFile))
  {
    std::cerr << "Input file path does not have a file name or the file does not exist." << std::endl;
    return -1;
  }

  if(!outFile.has_filename())
  {
    std::cerr << "Output file path does not have a file name." << std::endl;
    return -2;
  }

  if(offset < 0 || offset > 3)
  {
    std::cerr << "Invalid reading frame value." << std::endl;
    return -3;
  }

  ////////// opening files
  std::ifstream input(inFile);
  std::ofstream output(outFile, std::ios::trunc);

  if(!input.is_open() || !output.is_open())
  {
    std::cerr << "Error while openening input and/or output file." << std::endl;
    return -4;
  }

  ////////// tune to the requested coding frame
  if(offset > 0)
  {
    input.ignore(offset);
  }

  ////////// decode
  return decode(input, output, waitForStart, useEndCodon);
}