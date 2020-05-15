# Rosalind: Translating RNA into Protein
## Problem statement
### Problem
The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

*Given*: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

*Return*: The protein string encoded by s.

### Sample Dataset
AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
### Sample Output
MAMAPRTEINSTRING

## Requrement analysis
### Input data issues

#### The first codon
In the example, the first codon is the start codon (AUG), but it is not specified if all possible input data begin in the same way.
##### Possible solutions
1. Decode everything
  - pros: It's simple and it is not specified if the decoder should wait for the start codon, so maybe it is not needed.
  - cons: This is not how biological decoder works
2. Ignore codons before the start codon
  - pros: This is how the biological decoder works
  - cons: It requires more implementation
##### Implemented approach
Both (1) and (2) are implemented and (2) is the default behaviour. The behaviour can be changed using '--waitforstart' parameter.

#### The stop codons
In the example, the last codon is a stop codon (UGA), but it is not specified if all possible input data ends in the same way.
##### Possible solutions
1. Decode everything
  - pros: It's simple and it is not specified if the decoder should end after an en end codon, so maybe it is not needed.
  - cons: This is not how biological decoder works
2. Ignore codons after an end codon
  - pros: This is how the biological decoder works
  - cons: It requires more implementation
##### Implemented approach
Both (1) and (2) are implemented and (2) is the default behaviour. The behaviour can be changed using '--usestop' parameter. If "usestop" is true and the input data does not contain an end codon after the last start codon, the string is still decoded, but an error message is displayed.

#### Multiple start and stops
In the example, there is only one gene coded, but in the real sequences, multiple genes are coded in the same string. How the program should behave in that case?
##### Possible solutions
1. Decode everything
  - pros: It's simple to run and decode everything
  - cons: This is not how biological decoder works
2. Decode only between the start and stop codons. In the case of "start start stop" use the first start.
  - pros: This is how the biological decoder works
  - cons: It requires more implementation
3. Ignore everything after the first stop codon
  - pros: It's simple to run and decode everything
  - cons: This is not how biological decoder works
##### Implemented approach
Both (1) and (2) are implemented and (2) is the default behaviour. The behaviour can be changed using '--usestop' parameter.

#### RNA reading frames
A reading frame is defined by the initial triplet of nucleotides from which translation starts. It sets the frame for a run of successive, non-overlapping codons, which is known as an "open reading frame" (ORF). For example, the string 5'-AAATGAACG-3' (see https://en.wikipedia.org/wiki/Genetic_code#/media/File:Homo_sapiens-mtDNA~NC_012920-ATP8+ATP6_Overlap.svg), if read from the first position, contains the codons AAA, TGA, and ACG; if read from the second position, it contains the codons AAT and GAA; and if read from the third position, it contains the codons ATG and AAC. Every sequence can, thus, be read in its 5' → 3' direction in three reading frames, each producing a possibly distinct amino acid sequence. Source: https://en.wikipedia.org/wiki/Genetic_code
##### Possible solutions
1. Use only one frame
  - pros: It's simple
  - cons: This is not how biological decoder works
2. Generate multiple sequences for all coding frames
  - pros: This approach gives all possible outputs
  - cons: The example output does not contain multiple lines
3. Define the frame as an input parameter
  - pros: It gives a possibility to decode using all possible frames and outputs a single option
  - cons: It requires more implementation
##### Implemented approach
The (3) solution is implemented and the "O" frame is the default one. The behaviour can be changed using '--frame' parameter.

#### Invalid characters
According to the task description, it cannot happen that the input file contains invalid (other than A, U, C, G) character.
##### Implemented approach
If it happens anyway, the program immediately returns an error code and displays an error message. It may happen that in the output file are partially decoded data.

## Algorithm

### Description
The implemented solution reads the input file codon after codon and wites the decoded protein symbol immediately to the output file.
### Time complexity
The approach has to read all the data, so it has O(N) complexity, where N is the input data length.
### Space complexity
The approach has constant space complexity because only one codon is the buffer.

## Compilation
The code can be compiled using g++ 8.0 or newer.
`g++ -o rosalind_prot  main.cpp -lstdc++fs`

### Help
      ` Options: 
      -h | --help           Print this help 
      -i | --input          Input file (default in.txt) 
      -o | --output         Output file (default out.txt) 
      -f | --frame          RNA reading frame. See: https://en.wikipedia.org/wiki/Genetic_code  
                            (chapter: Features, subchapter: Reading frame) Valid: 0, 1, 2 (default: 0) 
      -w | --waitforstart   Flag if proteins decoded before the first start codon should be ignored 
                            Valid: 0, 1 (default: 1) << std:: endl <<
      -s | --usestop        Flag if proteins decoded after stop codon and before the next start codon  
                            should be ignored. Valid: 0, 1 (default: 1)  
      
      The application reads the RNA sequence from the input file and writes the decoded proteins to the output file 
      Task description: http://rosalind.info/problems/prot/ 
      The program returns 0 on success and an error code otherwise. 
      
      Error code | Meaning  
      -1         | Input file path does not have a filename or the file does not exist. 
      -2         | Output file path does not have a filename. 
      -3         | Invalid reading frame value 
      -4         | Error while openening input and/or output file. 
      -5         | Input file contains an invalid codon - a codon that contains characters other than A,U,C,G 
      -6         | Missing end codon - the decoded sequence in the output file is not finished 
      -7         | Error on parsing input parameters << std::endl`
      
### Test cases
Folder "tests" contains example inputs and outputs of the program.
