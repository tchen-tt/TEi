#ifndef __TESEQ__
#define __TESEQ__
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <Rcpp.h>
#include <string>
#include <cstring>
#include <map>
#include <vector>

using namespace Rcpp;
using std::map;
using std::strtok;


std::string maps = "NACNGNNNTNNNNNNN";
const char *int2char = maps.data();

map<char, char> transform = { {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'} };

#define SOFTCLIP_MATCH 1
#define MATCH_SOFTCLIP 2

typedef std::vector<char *> char_vector;

//  Process bam file.
void processSequence(bam1_t *b, bam_hdr_t *hdr, FILE *f, int length);
void processbam(std::string bamfile, double quantile, int length);
void processSequence1(bam1_t *b, bam_hdr_t *hdr, FILE *f, int length, int tsd);



// Convernt bam file to bed files.
int sortFunction(const void *x, const void *y);
int countNumber(int32_t *tids, int l);
int countTable(int32_t *tids, int pos, int l);
void bam2bed(bam1_t **read, int l, FILE *f, double threshold, bam_hdr_t *heaader);
void processSam2bed(std::string alignmentfile, std::string outbedfile, float ratio, int alignment);

#endif
