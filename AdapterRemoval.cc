/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 * Copyright (C) 2011 by Stinus Lindgreen                                *
 * stinus@binf.ku.dk                                                     *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
\*************************************************************************/

#include <math.h>
#include <sstream>
#include <cctype>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include <ctime>

using namespace std;

string VERSION="AdapterRemoval ver. 1.5.4";
string HELPTEXT="This program searches for and removes remnant adapter sequences from your read data. The program can analyze both single end and paired end data. Usage:\nAdapterRemoval [--file|file1 filename] [--file2 filename] [--basename filename] [--trimns] [--maxns max] [--trimqualities] [--minquality minimum] [--collapse] [--stats] [--version] [--mm mismatchrate] [--minlength len] [--minalignmentlength len] [--qualitybase base] [--shift num] [--pcr1 sequence] [--pcr2 sequence] [--5prime sequence] [--output1 file] [--output2 file] [--outputstats file] [--singleton file] [--singletonstats file] [--outputcollapsed file] [--outputcollapsedtruncated file] [--discarded file] [--settings file]\n\nFor detailed explanation of the parameters, please refer to the man page.\nIf nothing else, at least input your read data to the program (either stdin or from af fastq file).\nFor comments, suggestions and feedback please contact Stinus Lindgreen (stinus@binf.ku.dk).\n\nIf you use the program, please cite the paper:\nS. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337\nhttp://www.biomedcentral.com/1756-0500/5/337/\n";


string FIVEPRIME="";
string PCR1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
string PCR2="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
string shortmultiplex="AGACGTGTGCTCTTCCGATCT";
string genericindex="CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTC";

bool SCORINGMATRIX=false; // Use scoring matrix made using quality scores
int minimumlength=15; // The minimum genomic length we're interested in
int minalignmentlength=11; // The minimum required genomic overlap before collapsing reads into one

double matches[41][41];
double mismatches[41][41];
double globalmm=-1.0;
double singlemm=1.0/3.0;
double pairmm=0.1;

int qualitybase=33;
bool trimquals=false;
int minquality=2;
string minqual="";

bool removeNs=false;
int maxNsallowed=1000;
bool DEBUG=false;
bool collapse=false;
int shift=2; // Allow for slipping basepairs

int sumofqualities=0;

unsigned long int numberoffulllengthcollapsed=0;
unsigned long int numberoftruncatedcollapsed=0;
unsigned long int numberofseqswithadapter=0;
unsigned long int totalnumberofnucleotides=0;
unsigned long int totalnumberofgoodreads=0;

int noalignment=0;
int alignedok=0;
int alignedcrappy=0;
int keep1=0;
int discard1=0;
int keep2=0;
int discard2=0;

// Do we have two files?
bool usefile2=false;

// Should we write stats to separate output files?
bool printstats=false;

double errormatrix[51][2];

void calc_error_matrix(){
  double Perror,Ptrue;
  for(int i=0;i<=50;i++){
    Perror=pow(10,-0.1*(i+1.0));
    Ptrue=1.0-Perror;
    // Transform to log space
    errormatrix[i][0]=log(Ptrue);
    errormatrix[i][1]=log(Perror/3.0);
  }
}

string itoa(int a){
  string res="";
  int remainder;
  if(a == 0) res="0";
  while(a>0){
    remainder=a%10;
    a=int(a/10);
    res=char(remainder+48)+res;
  }
  return res;
}

void calc_matrices(){
  double PT1,PT2,PE1,PE2,PMATCH,PMM;
  for(int i=0;i<41;i++){
    if(i==40) PE1=0;
    else PE1=pow(10,-0.1*(i+1.0));
    PT1=1.0-PE1;
    //    cout<<"P_error("<<char(i+qualitybase+1)<<")="<<PE1<<"\n";
    for(int j=0;j<41;j++){
      if(j==40) PE2=0;
      else PE2=pow(10,-0.1*(j+1.0));
      PT2=1.0-PE2;
      PMATCH=PT1*PT2+(PE1*PE2)/3.0;
      PMM=1.0-PMATCH;
      matches[i][j]=log(PMATCH/0.25)/log(2);
      mismatches[i][j]=log(PMM/0.25)/log(2);
    }
  }
}

size_t min(size_t a,size_t b){
  if(a < b) return a;
  else return b;
}

inline string toupper(const string& s){
  string news = s;
  for(size_t i=0;i<s.length();i++){
    if(s[i] >= 'a' && s[i] <= 'z'){
      news[i] -= 32;
  }
}

  return news;
  }

inline string reverse(const string& seq){
  return string(seq.rbegin(), seq.rend());
}

string complement(const string& seq){
  string complement = seq;
  for(size_t i=0;i<seq.length();i++){
    if(seq[i]=='A')
      complement[i]='T';
    else if(seq[i]=='C')
      complement[i]='G';
    else if(seq[i]=='G')
      complement[i]='C';
    else if(seq[i]=='T')
      complement[i]='A';
  }
  return complement;
}

inline string replace_dots(const string& s){
  string news = s;
  for(size_t i=0;i<news.length();i++){
    if(news.at(i) == '.'){
      news.at(i) = 'N';
    }
  }

  return news;
}

double score(char n1,char n2,char q1, char q2){
  if(SCORINGMATRIX){
    // Use the scoring matrix built on the quality scores
    if(n1 == n2 && n1 != 'N') return matches[int(q1)-34][int(q2)-34];
    else return mismatches[int(q1)-34][int(q2)-34];
  }
  else{
    if(n1 == 'N' || n2 == 'N') return 0;
    else if(n1 == n2) return 1;
    else return -1;
  }
}

string reverse_complement(const string& seq){
  return reverse(complement(seq));
}

int trimNs(string& sequence,string& qualities,string& stats){
  int length=0;
  int begin=sequence.find_first_not_of('N');
  int trimmed3=0;
  int trimmed5=0;
  if(DEBUG) cerr<<"Trimming Ns\nInitial start: "<<begin<<"\n";
  if(begin>-1){
    length=(sequence.find_last_not_of('N'))-begin+1;
    trimmed5=begin;
    trimmed3=sequence.length()-length-trimmed5;
  }
  else{
    // All Ns
    begin=sequence.length();
    length=0;
    trimmed5=0;
    trimmed3=sequence.length();
  }
  if(DEBUG){
    cerr<<"Final start: "<<begin<<"\tLength: "<<length<<"\n";
    cerr<<"Found Ns on first "<<trimmed5<<" positions and on last "<<trimmed3<<" positions in:\n"<<sequence<<"\n"<<qualities<<"\n";
  }
  stats=stats+"\t5pN:"+itoa(trimmed5)+"\t3pN:"+itoa(trimmed3); 
  sequence=sequence.substr(begin,length);
  qualities=qualities.substr(begin,length);
  if(DEBUG){
    cerr<<"After trimming: \n"<<sequence<<"\n"<<qualities<<"\n";
  }
  return trimmed3+trimmed5;
}

int trimqualities(string& sequence, string& qualities, string& stats){
  // Trim reads by removing low qualities as specified in string minqual.
  // Returns number of positions trimmed in 3' end
  int length=0;
  int begin=qualities.find_first_not_of(minqual);
  int trimmed3=0;
  int trimmed5=0;
  if(DEBUG) cerr<<"Trimming qualities\nInitial start: "<<begin<<"\n";
  if(begin>-1){
    length=(qualities.find_last_not_of(minqual))-begin+1;
    trimmed5=begin;
    trimmed3=qualities.length()-length-trimmed5;
  }
  else{
    // Whole sequence is of low quality
    begin=qualities.length();
    length=0;
    trimmed5=0;
    trimmed3=qualities.length();
  }
  if(DEBUG){
    cerr<<"Final start: "<<begin<<"\tLength: "<<length<<"\n";
    cerr<<"Found "<<minqual<<" on first "<<trimmed5<<" positions and on last "<<trimmed3<<" positions in:\n"<<sequence<<"\n"<<qualities<<"\n";
  }
  stats=stats+"\t5pQ:"+itoa(trimmed5)+"\t3pQ:"+itoa(trimmed3);
  sequence=sequence.substr(begin,length);
  qualities=qualities.substr(begin,length);
  if(DEBUG){
    cerr<<"After trimming: \n"<<sequence<<"\n"<<qualities<<"\n";
  }
  return trimmed3+trimmed5;
}

void collapsealignment(string& seq1,string& seq2,string& qual1, string& qual2){
  if(DEBUG) cerr<<"Collapse: \n"<<seq1<<"\t"<<qual1<<"\n"<<seq2<<"\t"<<qual2<<"\n"<<"Qual1\tPtrue1\tPerror1\tQual2\tPtrue2\tPerror2\n";

  char nucleotides[4]={'A','C','G','T'};
  double maxprob=0.9999;
  double Ptrue,Perror;
  vector<vector<double> > probabilities;
  
  // Calculate probability matrix based on the two reads and their qualities
  for(size_t i=0;i<seq1.length();i++){
    vector<double> probs(4,0.0); 
    if(seq1[i] != 'N'){
      Ptrue=errormatrix[int(qual1[i])-qualitybase][0];
      Perror=errormatrix[int(qual1[i])-qualitybase][1];
      for(int j=0;j<4;j++){
	if(seq1[i]==nucleotides[j]) probs[j]=Ptrue;
	else probs[j]=Perror;
      }
    }
    if(DEBUG) cerr<<qual1[i]<<"\t"<<Ptrue<<"\t"<<Perror<<"\t";
    probabilities.push_back(probs);
    if(seq2[i] != 'N'){
      Ptrue=errormatrix[int(qual2[i])-qualitybase][0];
      Perror=errormatrix[int(qual2[i])-qualitybase][1];
      for(int j=0;j<4;j++){
	if(seq2[i]==nucleotides[j]) probabilities[i][j] += Ptrue;
	else  probabilities[i][j] += Perror;
      }
    }
    if(DEBUG) cerr<<qual2[i]<<"\t"<<Ptrue<<"\t"<<Perror<<"\n";
  }

  // Find new consensus sequence
  for(size_t i=0;i<probabilities.size();i++){
    double largest=probabilities[i][0];
    char nuc='A';

    // Handle the case of different nucleotides with identical quality scores
    if(seq1[i] != seq2[i] && qual1[i] == qual2[i]){
      // Solution 1: Choose the first nucleotide...
      // Solution 2: Use an N and set low quality
      // Solution 3: Choose one of the nucleotides at random (if one is an N: Always use the other)
      if(seq1[i] == 'N' || seq2[i] == 'N'){
	seq1[i]=(seq1[i]=='N')?seq2[i]:seq1[i];
      }
      else{
	int shuffle=rand()%2;
	seq1[i]=(shuffle)?seq1[i]:seq2[i];
      }
    }
    else{
      for(int j=1;j<4;j++){
	if(probabilities[i][j]>largest){
	  largest=probabilities[i][j];
	  nuc=nucleotides[j];
	}
      }
      
      double normconstant=0.0;
      for(int j=0;j<4;j++) normconstant += exp(probabilities[i][j]-largest);
      normconstant=largest+log(normconstant);
      double normprob=exp(largest-normconstant);
      if( (normprob-0.25)<0.0000001 && (normprob-0.25)>-0.0000001 ) nuc='N';
      seq1[i]=nuc;
      if( normprob >= maxprob)
	qual1[i]=char(qualitybase+41);
      else
	qual1[i]=char(-10*log10(1.0-normprob)+qualitybase-0.5);
      //	qual1[i]=char(-10*log10(1.0-normprob)+qualitybase+0.5);
    }
  }
  if(DEBUG) cerr<<"Collapse: \n"<<seq1<<"\t"<<qual1<<"\n";
}


void trim5primeend(string& read,string& qualities){
  // Quick and dirty. Allow alignment to be shifted by up to <shift> positions. Accept 1 mismatch.
  int mmthreshold=1;
  int bestmm=mmthreshold+1;
  int bestshift=-1;
  int trimmed5=0;
  for(int s=0;s<=shift;s++){
    int mismatches=0;
    for(size_t i=0;i<FIVEPRIME.length()-s;i++){
      if(FIVEPRIME[i+s]!=read[i]) mismatches++;
    }
    if(mismatches<bestmm){
      bestmm=mismatches;
      bestshift=s;
    }
  }
  if(bestmm<=mmthreshold){
    read=read.substr(FIVEPRIME.length()-bestshift);
    qualities=qualities.substr(FIVEPRIME.length()-bestshift);
    trimmed5=FIVEPRIME.length()-bestshift;
  }
}

// Performs the alignment of reads and primers. Reads and qualities are truncated and various statistics are calculated.
void pairwisealignment(string read1,string read2,string qualities1,string qualities2,string id1,string id2,ostream& output1,ostream& output2,ostream& discarded, ostream& singleton, ostream& collapsed, ostream& collapsedtruncated, ostream& singletonstats, ostream& truncatedstats){
  // Align read1 to reverse-complement read2
  int trimmed=0;
  // Concatenate PCR adapters and reads
  int diff;
  string seq1,seq2,qual1,qual2,stats1="",stats2="";
  if(usefile2){
    seq1=PCR2+read1;
    qual1=string(PCR2.length(),char(qualitybase+41))+qualities1; 
    diff=(int)(read1.length()-PCR2.length());
    if(diff>0){
      seq1=string(diff,'N')+seq1;
      qual1=string(diff,char(qualitybase+41))+qual1;
    }
    seq2=reverse_complement(read2)+PCR1;
    qual2=reverse(qualities2)+string(PCR1.length(),char(qualitybase+41));
    // Pad the sequencese with Ns to get equal length
    diff=(int)(read2.length()-PCR1.length());
    if(diff>0){
      seq2=seq2+string(diff,'N');
      qual2=qual2+string(diff,char(qualitybase+41));
    }  
    diff=(int)(seq2.length()-seq1.length());
    if(diff>0){
      seq1=string(diff,'N')+seq1;
      qual1=string(diff,char(qualitybase+41))+qual1;
    }
    else{
      seq2=seq2+string(-1*diff,'N');
      qual2=qual2+string(-1*diff,char(qualitybase+41));
    }
  }
  else{
    seq1=read1;
    qual1=qualities1;

    diff=read1.length()+shift-read2.length();
    if(diff<0){
      // Single end case, adapter longer than read: Truncate adapter
      seq2=read2.substr(0,read1.length()+shift);
      qual2=qualities2.substr(0,read1.length()+shift);
    }
    else{
      // Single end case, read longer than adapter: Extend with Ns
      seq2=read2+string(diff,'N');
      qual2=qualities2+string(diff,char(qualitybase+41));
    }
  }

  // Initialize the matrices
  if(DEBUG && seq1.length() != seq2.length()) cerr<<"Read1 length: "<<seq1.length()<<" Read2 length: "<<seq2.length()<<"\n";

  double matrix[seq1.length()+1][seq2.length()+1];
  int mmmatrix[seq1.length()+1][seq2.length()+1];
  int nmatrix[seq1.length()+1][seq2.length()+1];

  for(size_t i=0;i<=seq1.length();i++){
    matrix[i][0]=0.0;
    mmmatrix[i][0]=0;
    nmatrix[i][0]=0;
  }
  for(int j=0;j<=shift;j++){
    matrix[0][j]=0;
    mmmatrix[0][j]=0;
    nmatrix[0][j]=0;
  }
  
  // Calculate upper matrix
  for(size_t j=1;j<=seq2.length();j++){
    for(size_t i=max(1,(int)j-shift);i<=seq1.length();i++){
      if(seq1[i-1] == 'N' || seq2[j-1]=='N'){
	// Align to N - neutral
	matrix[i][j]=matrix[i-1][j-1];
	mmmatrix[i][j]=mmmatrix[i-1][j-1];
	nmatrix[i][j]=nmatrix[i-1][j-1]+1;
      }
      else if(seq1[i-1] == seq2[j-1]){
	// Match - score+1
	matrix[i][j]=matrix[i-1][j-1]+1;	  
	mmmatrix[i][j]=mmmatrix[i-1][j-1];
	nmatrix[i][j]=nmatrix[i-1][j-1];
      }
      else{
	// Mismatch - score-1
	matrix[i][j]=matrix[i-1][j-1]-1;
	mmmatrix[i][j]=mmmatrix[i-1][j-1]+1;
	nmatrix[i][j]=nmatrix[i-1][j-1];
      }
    }
  }

  if(DEBUG){
    cerr<<seq1<<"\n"<<qual1<<"\n"<<seq2<<"\n"<<qual2<<"\nScore matrix:\n";
    for(size_t j=0;j<=seq2.length();j++){
      for(size_t i=0;i<=seq1.length();i++){
	if(i>=j-shift)
	  cerr<<matrix[i][j]<<"\t";
	else
	  cerr<<"-\t";
      }
      cerr<<"\n";
    }
    cerr<<"\nMismatch matrix:\n";
    for(size_t j=0;j<=seq2.length();j++){
      for(size_t i=0;i<=seq1.length();i++){
	if(i>=j-shift)
	  cerr<<mmmatrix[i][j]<<"\t";
	else
	  cerr<<"-\t";
      }
      cerr<<"\n";
    }
    cerr<<"\nN matrix\n";
    for(size_t j=0;j<=seq2.length();j++){
      for(size_t i=0;i<=seq1.length();i++){
	if(i>=j-shift)
	  cerr<<nmatrix[i][j]<<"\t";
	else
	  cerr<<"-\t";
      }
      cerr<<"\n";
    }
    cerr<<"\n";
  }

  // Find best match
  // By looking in the rightmost column - i.e. for i=seq1.length() - and using 1<=j<=seq2.length() 
  // we examine all alignments between the 3' end of seq1 and the 5' end of seq2.
  // Using matrix[seq1.length][1] corresponds to aligning 1 nucleotide from each sequence.
  // Using matrix[seq1.length][seq2.length] corresponds to aligning the two in full length,
  // and this also means we have no genomic DNA.
  // This is the longest, highest scoring alignment with an acceptable number of mismatches. 

  int alignmentlength=0;
  int maxScore=0;
  int num_of_mm=0;
  bool unaligned=true;
  int debugj=-1;
  int mmthreshold;
  int numberofNs=0,numberofNs1=0,numberofNs2=0;
  for(size_t j=1;j<=seq2.length();j++){
    if(j<6) mmthreshold=0;
    else if(j<10) mmthreshold=1;
    else mmthreshold=(int)(globalmm*(double)(j-nmatrix[seq1.length()][j]));
    //if(matrix[seq1.length()][j] >= maxScore && mmmatrix[seq1.length()][j]<=mmthreshold){
    if(matrix[seq1.length()][j] >= maxScore){
      unaligned=false;
      alignmentlength=j;
      maxScore=matrix[seq1.length()][j];
      num_of_mm=mmmatrix[seq1.length()][j];
      numberofNs=nmatrix[seq1.length()][j];
      debugj=j;
    }
  }
  //  if(usefile2 || alignmentlength>10) mmthreshold=(int)(globalmm*(double)alignmentlength);
  if( alignmentlength < 6 ) mmthreshold=0;
  else if(alignmentlength < 10) mmthreshold=1;
  else mmthreshold=(int)(globalmm*(double)(alignmentlength-numberofNs));
  // If too many mismatches
  if(num_of_mm > mmthreshold) unaligned=true;

  // If the aligned part is too short to collapse the reads, treat them as unaligned
  if(collapse && (alignmentlength-numberofNs)<minalignmentlength) unaligned=true;

  if(DEBUG) cerr<<id1<<"\tUnaligned? "<<unaligned<<"\tLength: "<<alignmentlength<<"\tEffective length: "<<(alignmentlength-numberofNs)<<"\tNumber of Ns: "<<numberofNs<<"\tMM: "<<num_of_mm<<"\tMM threshold: "<<mmthreshold<<"\tScore: "<<maxScore<<"\tj: "<<debugj<<"\n";

  if(unaligned){
    // Could not align the two reads - output unaltered reads
    if(removeNs){
      trimNs(read1,qualities1,stats1);
      if(usefile2) trimNs(read2,qualities2,stats2);
    }
    if(trimquals){
      trimqualities(read1,qualities1,stats1);
      if(usefile2) trimqualities(read2,qualities2,stats2);
    }

    // Count number of Ns
    numberofNs1=0;
    numberofNs2=0;
    for(size_t i=0;i<=read1.length();i++) if(read1[i]=='N') numberofNs1++;
    if(usefile2) for(size_t i=0;i<=read2.length();i++) if(read2[i]=='N') numberofNs2++;

    // Check length and number of Ns
    if( (read1.length() >= (size_t) minimumlength && numberofNs1 <= maxNsallowed) && ( read2.length() >= (size_t) minimumlength && numberofNs2 <= maxNsallowed) ){
      noalignment++;
      if(printstats) truncatedstats<<id1<<stats1<<"\tUnaligned";
      output1<<id1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
      totalnumberofnucleotides += read1.length();
      totalnumberofgoodreads++;
      if(usefile2){
	if(printstats) truncatedstats<<"\t"<<id2<<stats2<<"\tUnaligned";
	output2<<id2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	totalnumberofnucleotides += read2.length();
	totalnumberofgoodreads++;
      }
      if(printstats) truncatedstats<<"\n";
    }
    else if( (read1.length() >= (size_t) minimumlength && numberofNs1 <= maxNsallowed) && (read2.length() < (size_t) minimumlength ||  numberofNs2 > maxNsallowed) ){
      keep1++;
      if(printstats) singletonstats<<id1<<stats1<<"\tUnaligned\n";
      totalnumberofnucleotides += read1.length();
      totalnumberofgoodreads++;
      // Keep read1. If we use only 1 file, output to output1. Otherwise, output to singleton
      if(usefile2){
	discard2++;
	singleton<<id1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
	replace(stats2.begin(),stats2.end(),'\t','_');
	discarded<<id2<<"_"<<stats2<<"_unaligned_tooshort_Ns:"<<numberofNs2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
      }
      else{
	output1<<id1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
      }
    }
    else if( (read1.length() < (size_t) minimumlength || numberofNs1 > maxNsallowed) && (read2.length() >= (size_t) minimumlength && numberofNs2 <= maxNsallowed) ){
      discard1++;
      replace(stats1.begin(),stats1.end(),'\t','_');
      discarded<<id1<<"_"<<stats1<<"_unaligned_tooshort_Ns:"<<numberofNs1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
      if(usefile2){
	keep2++;
	singleton<<id2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	totalnumberofnucleotides += read2.length();
	totalnumberofgoodreads++;
	if(printstats) singletonstats<<id2<<stats2<<"\tUnaligned\n";
      }
    }
    else{ // Both too short
      discard1++;
      replace(stats1.begin(),stats1.end(),'\t','_');
      discarded<<id1<<"_"<<stats1<<"_unaligned_tooshort_Ns:"<<numberofNs1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
      if(usefile2){
	discard2++;
	replace(stats2.begin(),stats2.end(),'\t','_');
	discarded<<id2<<"_"<<stats2<<"_unaligned_tooshort_Ns:"<<numberofNs2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
      }
    }
  }
  else{ // Aligned
    int genomiclength=0;
    int adapterlength=0;
    if(usefile2){
      if((size_t) alignmentlength <= min(read1.length(),read2.length())){
	// We have only genomic overlap and no adapter
	genomiclength=read1.length()+read2.length()-alignmentlength;
	// The is the length of the combined read
      }
      else{
	//We have some adapter sequence in the reads
	adapterlength=alignmentlength-read1.length();
	genomiclength=alignmentlength-2*adapterlength;
      }
    }
    else{
      adapterlength=alignmentlength;
      genomiclength=max(0,(int)seq1.length()-alignmentlength);
    }
     
    if(DEBUG){
      cerr<<"Alignmentlength: "<<alignmentlength<<" Adapter length: "<<adapterlength<<" Genomic length: "<<genomiclength<<"\n";
    }

    // Now, create new truncated reads based on adapters
    string newseq1,newseq2,newqual1,newqual2;
    if(!collapse){
      stats1="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm);
      newseq1=read1.substr(0,max(0,(int)read1.length()-adapterlength));
      newqual1=qualities1.substr(0,max(0,(int)read1.length()-adapterlength));
      if(usefile2){
	stats2="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm);
	newseq2=read2.substr(0,max(0,(int)read2.length()-adapterlength));
	newqual2=qualities2.substr(0,max(0,(int)read2.length()-adapterlength));
      }
      else{
	newseq2=read2;
	newqual2=qualities2;
      }
    }
    else{
      // Remember to work on original read1 but reverse complement of read2!
      if(alignmentlength < genomiclength){
	// We combine the two reads with a new consensus in the middle
	string overlapseq1=read1.substr(read1.length()-alignmentlength);
	string overlapseq2=seq2.substr(0,alignmentlength);
	string overlapqual1=qualities1.substr(qualities1.length()-alignmentlength);
	//	string overlapqual2=qualities2.substr(0,alignmentlength);
	string overlapqual2=qual2.substr(0,alignmentlength);
	if(DEBUG) cerr<<"Before collapse 1: Alignmentlength: "<<alignmentlength<<" Read1 length: "<<read1.length()<<" Read2 length: "<<read2.length()<<"\n"<<overlapseq1<<"\n"<<overlapqual1<<"\n"<<overlapseq2<<"\n"<<overlapqual2<<"\n";
	collapsealignment(overlapseq1,overlapseq2,overlapqual1,overlapqual2);
	if(DEBUG) cerr<<"After collapse 1: Alignmentlength: "<<alignmentlength<<" Read1 length: "<<read1.length()<<" Read2 length: "<<read2.length()<<"\n"<<overlapseq1<<"\n"<<overlapqual1<<"\n"<<overlapseq2<<"\n"<<overlapqual2<<"\n";
	stats1="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm)+"\tCollapsed";
	stats2="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm)+"\tCollapsed";
	newseq1=read1.substr(0,read1.length()-alignmentlength)+overlapseq1+seq2.substr(alignmentlength,read2.length()-alignmentlength);
	newqual1=qualities1.substr(0,read1.length()-alignmentlength)+overlapqual1+qual2.substr(alignmentlength,read2.length()-alignmentlength);
	if(DEBUG) cerr<<"Sequence length: "<<newseq1.length()<<" Quality length: "<<newqual1.length()<<"\n"<<id1<<"\n"<<newseq1<<"\n"<<newqual1<<"\n";
      }
      else{
	// We cut off the adapter and combine the two reads into one shorter consensus
	stats1="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm)+"\tCollapsed";
	newseq1=read1.substr(0,genomiclength);
	newqual1=qualities1.substr(0,genomiclength);
	stats2="\talignment:"+itoa(alignmentlength)+"\tadapter:"+itoa(adapterlength)+"\tMM:"+itoa(num_of_mm)+"\tCollapsed";
	newseq2=seq2.substr(adapterlength,genomiclength);
	newqual2=qual2.substr(adapterlength,genomiclength);
	if(DEBUG) cerr<<"Before collapse 2:\n"<<newseq1<<"\n"<<newqual1<<"\n"<<newseq2<<"\n"<<newqual2<<"\n";
	collapsealignment(newseq1,newseq2,newqual1,newqual2);
	if(DEBUG) cerr<<"After collapse 2:\n"<<newseq1<<"\n"<<newqual1<<"\n"<<newseq2<<"\n"<<newqual2<<"\n";
      }
    }
    // Trim Ns and qualities if specified
    if(removeNs){
      trimmed += trimNs(newseq1,newqual1,stats1);
      if(usefile2 && !collapse) trimmed += trimNs(newseq2,newqual2,stats2);
    }
    if(trimquals){
      trimmed += trimqualities(newseq1,newqual1,stats1);
      if(usefile2 && !collapse) trimmed += trimqualities(newseq2,newqual2,stats2);
    }

    // count number of Ns in trimmed reads
    numberofNs1=0;
    numberofNs2=0;
    for(size_t i=0;i<=newseq1.length();i++) if(newseq1[i]=='N') numberofNs1++;
    if(usefile2) for(size_t i=0;i<=newseq2.length();i++) if(newseq2[i]=='N') numberofNs2++;

    //The two reads were aligned - but is the quality good enough?
    //    if( ( newseq1.length() >= minimumlength && numberofNs1 <= maxNsallowed ) && ( newseq2.length() >= minimumlength && numberofNs2 <= maxNsallowed) ){
    if( ( newseq1.length() >= (size_t) minimumlength && numberofNs1 <= maxNsallowed ) && (collapse || ( newseq2.length() >= (size_t) minimumlength && numberofNs2 <= maxNsallowed) ) ){
      // Both reads are of adequate length
      if(num_of_mm <= mmthreshold && maxScore > 0){
	// Lengths ok and the score plus number of mismatches are both ok
	alignedok++;
	if(collapse){
	  if(adapterlength>0) numberofseqswithadapter++;
	  if(trimmed==0){
	    collapsed<<"@M_"<<id1.substr(1)<<"\n"<<newseq1<<"\n+\n"<<newqual1<<"\n";
	    numberoffulllengthcollapsed++;
	  }
	  else{
	    // Collapsed and trimmed in 5' and/or 3' end
	    collapsedtruncated<<"@MT_"<<id1.substr(1)<<"\n"<<newseq1<<"\n+\n"<<newqual1<<"\n";
	    numberoftruncatedcollapsed++;
	  }
	  totalnumberofnucleotides += newseq1.length();
	  totalnumberofgoodreads++;
	  if(printstats) singletonstats<<id1<<"\t"<<id2<<stats1<<"\n";
	}
	else{
	  // Not collapse
	  output1<<id1<<"\n"<<newseq1<<"\n+\n"<<newqual1<<"\n";
	  totalnumberofnucleotides += newseq1.length();
	  totalnumberofgoodreads++;
	  if(adapterlength>0) numberofseqswithadapter++;
	  if(printstats){
	    truncatedstats<<id1<<stats1<<"\tAligned";
	    if(!usefile2) truncatedstats<<"\n";
	    else truncatedstats<<"\t"<<id2<<stats2<<"\tAligned"<<"\n";
	  }
	  if(usefile2){
	    output2<<id2<<"\n"<<newseq2<<"\n+\n"<<newqual2<<"\n";
	    totalnumberofnucleotides += newseq2.length();
	    totalnumberofgoodreads++;
	    if(adapterlength>0) numberofseqswithadapter++;
	  }
	}
      }
      else{
	// Lengths ok but the alignment is not good enough
	alignedcrappy++;
	stats1="";
	stats2="";
	if(removeNs){
	  trimmed += trimNs(read1,qualities1,stats1);
	  if(usefile2) trimmed += trimNs(read2,qualities2,stats2);
	}
	if(trimquals){
	  trimmed += trimqualities(read1,qualities1,stats1);
	  if(usefile2) trimmed += trimqualities(read2,qualities2,stats2);
	}

	// Count number of Ns
	numberofNs1=0;
	numberofNs2=0;
	for(size_t i=0;i<=read1.length();i++) if(read1[i]=='N') numberofNs1++;
	if(usefile2) for(size_t i=0;i<=read2.length();i++) if(read2[i]=='N') numberofNs2++;

	// Check length
	if( (read1.length() >= (size_t) minimumlength && numberofNs1 <= maxNsallowed) && (read2.length() >= (size_t) minimumlength && numberofNs2 <= maxNsallowed) ){
	  output1<<id1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
	  if(printstats) truncatedstats<<id1<<stats1<<"\tAlignmentFailed";
	  totalnumberofnucleotides += read1.length();
	  totalnumberofgoodreads++;
	  if(adapterlength>0) numberofseqswithadapter++;
	  if(usefile2){
	    output2<<id2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	    totalnumberofnucleotides += read2.length();
	    totalnumberofgoodreads++;
	    if(adapterlength>0) numberofseqswithadapter++;
	    if(printstats) truncatedstats<<"\t"<<id2<<stats2<<"\tAlignmentFailed";
	  }
	  if(printstats) truncatedstats<<"\n";
	}
	else if( (read1.length() >= (size_t) minimumlength && numberofNs1 <= maxNsallowed) && (read2.length() < (size_t) minimumlength || numberofNs2 > maxNsallowed) ){
	  keep1++;
	  singleton<<id1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
	  totalnumberofnucleotides += read1.length();
	  totalnumberofgoodreads++;
	  if(adapterlength>0) numberofseqswithadapter++;
	  if(printstats) singletonstats<<id1<<stats1<<"\tAlignmentFailed"<<"\n";
	  if(usefile2){
	    discard2++;
	    replace(stats2.begin(),stats2.end(),'\t','_');
	    discarded<<id2<<"_"<<stats2<<"_Ns:"<<numberofNs2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	  }
	}
	else if( (read1.length() < (size_t) minimumlength ||  numberofNs1 > maxNsallowed) && (read2.length() >= (size_t) minimumlength && numberofNs2 <= maxNsallowed) ){
	  discard1++;
	  replace(stats1.begin(),stats1.end(),'\t','_');
	  discarded<<id1<<"_"<<stats1<<"_Ns:"<<numberofNs1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
	  if(usefile2){
	    keep2++;
	    singleton<<id2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	    totalnumberofnucleotides += read2.length();
	    totalnumberofgoodreads++;
	    if(adapterlength>0) numberofseqswithadapter++;
	    if(printstats) singletonstats<<id2<<stats2<<"\tAlignmentFailed"<<"\n";
	  }
	}
	else{ // Both too short
	  discard1++;
	  replace(stats1.begin(),stats1.end(),'\t','_');
	  discarded<<id1<<"_"<<stats1<<"_AlignmentFailed_Ns:"<<numberofNs1<<"\n"<<read1<<"\n+\n"<<qualities1<<"\n";
	  if(usefile2){
	    discard2++;
	    replace(stats2.begin(),stats2.end(),'\t','_');
	    discarded<<id2<<"_"<<stats2<<"_AlignmentFailed_Ns:"<<numberofNs2<<"\n"<<read2<<"\n+\n"<<qualities2<<"\n";
	  }
	}
	if(DEBUG){	
	  cerr<<num_of_mm<<"\t"<<mmthreshold<<"\t"<<maxScore<<"\n";
	  cerr<<id1<<"\t"<<id1<<"\n"<<read1<<"\t"<<qualities1<<"\n"<<newseq1<<"\t"<<newqual1<<"\n";
	  if(usefile2) cerr<<id2<<"\t"<<id2<<"\n"<<read2<<"\t"<<qualities2<<"\n"<<newseq2<<"\t"<<newqual2<<"\n";
	}
      }
    }
    else{
      // One or both reads are too short
      if(newseq1.length() < (size_t) minimumlength || numberofNs1 > maxNsallowed){
	// Discard read1
	discard1++;
	replace(stats1.begin(),stats1.end(),'\t','_');
	discarded<<id1<<"_"<<stats1<<"_tooshort_Ns:"<<numberofNs1<<"\n"<<newseq1<<"\n+\n"<<newqual1<<"\n";
      }
      else{
	keep1++;
	singleton<<id1<<"\n"<<newseq1<<"\n+\n"<<newqual1<<"\n";
	totalnumberofnucleotides += newseq1.length();
	totalnumberofgoodreads++;
	if(adapterlength>0) numberofseqswithadapter++;
	if(printstats) singletonstats<<id1<<stats1<<"\n";
      }
      if(usefile2){
	if(newseq2.length() < (size_t) minimumlength || numberofNs2 > maxNsallowed){
	  // Discard read2
	  discard2++;
	  replace(stats2.begin(),stats2.end(),'\t','_');
	  discarded<<id2<<"_"<<stats2<<"_tooshort_Ns:"<<numberofNs2<<"\n"<<newseq2<<"\n+\n"<<newqual2<<"\n";
	}
	else{
	  keep2++;
	  singleton<<id2<<"\n"<<newseq2<<"\n+\n"<<newqual2<<"\n";
	  totalnumberofnucleotides += newseq2.length();
	  totalnumberofgoodreads++;
	  if(adapterlength>0) numberofseqswithadapter++;
	  if(printstats) singletonstats<<id2<<stats2<<"\n";
	}
      }
    }
  }
}

// Main function. Read in sequences. Remember names.
int main(int argc, char *argv[]){

  ios_base::sync_with_stdio(false);

  srand(time(NULL));
  calc_matrices();
  string adapter1,adapter2;

  bool usefile1=false;
  bool outfilegiven=false;
  char* file1;
  char* file2;
  string outfile="your_output";
  int lines=0;
  bool pcr1set=false;
  bool fiveprimeset=false;
  bool pcr2set=false;

  ostream *output1;
  ostream *output2;
  ostream *singleton;
  ostream *collapsed;
  ostream *collapsedtruncated;
  ostream *discarded;
  ostream *settings;
  ostream *truncatedstats;
  ostream *singletonstats;

  bool output1set=false;
  bool output2set=false;
  bool singletonset=false;
  bool collapsedset=false;
  bool collapsedtruncatedset=false;
  bool discardedset=false;
  bool settingsset=false;
  bool truncatedstatsset=false;
  bool singletonstatsset=false;

  for(int i=1;i<argc;i++){
    if(string(argv[i]) == "--pcr1") {PCR1=toupper(string(argv[i+1]));i++;pcr1set=true;}
    else if(string(argv[i]) == "--5prime") {FIVEPRIME=toupper(string(argv[i+1]));i++;fiveprimeset=true;}
    else if(string(argv[i]) == "--pcr2") {PCR2=toupper(string(argv[i+1]));i++;pcr2set=true;}
    else if(string(argv[i]) == "--basename") {outfile=string(argv[i+1]);outfilegiven=true;i++;}
    else if(string(argv[i]) == "--mm") {globalmm=atof(argv[i+1]);i++;}
    else if(string(argv[i]) == "--trimns") {removeNs=true;}
    else if(string(argv[i]) == "--maxns") {maxNsallowed=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--trimqualities") {trimquals=true;}
    else if(string(argv[i]) == "--qualitybase") {qualitybase=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--minquality") {minquality=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--minlength") {minimumlength=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--minalignmentlength") {minalignmentlength=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--file1" || string(argv[i]) == "--file") {file1=argv[i+1]; usefile1=true;i++;}
    else if(string(argv[i]) == "--file2") {file2=argv[i+1]; usefile2=true;i++;}
    else if(string(argv[i]) == "--debug") {DEBUG=true;}
    else if(string(argv[i]) == "--shift") {shift=atoi(argv[i+1]);i++;}
    else if(string(argv[i]) == "--collapse") {collapse=true;}
    else if(string(argv[i]) == "--stats") {printstats=true;}
    else if(string(argv[i]) == "--version") {cerr<<VERSION<<"\n";return 1;}

    else if(string(argv[i]) == "--discarded") {discarded=new ofstream(argv[i+1],ofstream::out); discardedset=true; i++;}
    else if(string(argv[i]) == "--settings") {settings=new ofstream(argv[i+1],ofstream::out); settingsset=true; i++;}
    else if(string(argv[i]) == "--output1") {output1=new ofstream(argv[i+1],ofstream::out); output1set=true; i++;}
    else if(string(argv[i]) == "--output2") {output2=new ofstream(argv[i+1],ofstream::out); output2set=true; i++;}
    else if(string(argv[i]) == "--singleton") {singleton=new ofstream(argv[i+1],ofstream::out); singletonset=true; i++;}
    else if(string(argv[i]) == "--outputcollapsed") {collapsed=new ofstream(argv[i+1],ofstream::out); collapsedset=true; i++;}
    else if(string(argv[i]) == "--outputcollapsedtruncated") {collapsedtruncated=new ofstream(argv[i+1],ofstream::out); collapsedtruncatedset=true; i++;}
    else if(string(argv[i]) == "--outputstats") {truncatedstats=new ofstream(argv[i+1],ofstream::out); truncatedstatsset=true; i++;}
    else if(string(argv[i]) == "--singletonstats") {singletonstats=new ofstream(argv[i+1],ofstream::out); singletonstatsset=true; i++;}
    else if(string(argv[i]) == "--help" || string(argv[i]) == "--h" || string(argv[i]) == "-help" || string(argv[i]) == "-h")  {cout<<VERSION<<"\n"<<HELPTEXT<<"\n"; return 0;}
    else {cerr<<VERSION<<"\n"<<"\n"<<"Unknown argument: "<<string(argv[i])<<"\n"<<"\n"<<HELPTEXT<<"\n"; return 1;}
  }
  // Check for invalid combinations of settings
  if(collapse && !usefile2){
    cerr<<"You tried to collapse a single file. Collapse has been reset to false.\n";
    collapse=false;
  }

  // Set the quality for trimming
  if(trimquals){
    if(minquality > qualitybase) // Assume we have been given the actual value, calculate offset
      minquality=minquality-qualitybase;
    for(int i=0;i<=minquality;i++)
      minqual=minqual+char(qualitybase+i);
  }

  // Create output files if needed
  // Always output sequences that fail to this file
  if(!discardedset) discarded=new ofstream((outfile+".discarded").c_str(),ofstream::out);
  // And always output program settings and overall statistics to this file
  if(!settingsset) settings=new ofstream((outfile+".settings").c_str(),ofstream::out);

  if(!usefile2){
    // Just one file - use stdout if basename not given
    if(!output1set){
      if(outfilegiven) output1=new ofstream((outfile+".truncated").c_str(),ofstream::out);
      else output1=&cout;
      if(printstats){
	if(!truncatedstatsset) truncatedstats=new ofstream((outfile+".truncated.stats").c_str(),ofstream::out);
      }
    }
  }
  else{
    // Two files - write to two output files
    if(!output1set) output1=new ofstream((outfile+".pair1.truncated").c_str(),ofstream::out);
    if(!output2set) output2=new ofstream((outfile+".pair2.truncated").c_str(),ofstream::out);
    if(!singletonset) singleton=new ofstream((outfile+".singleton.truncated").c_str(),ofstream::out);
    if(collapse){
      if(!collapsedset) collapsed=new ofstream((outfile+".collapsed").c_str(),ofstream::out);
      if(!collapsedtruncatedset) collapsedtruncated=new ofstream((outfile+".collapsed.truncated").c_str(),ofstream::out);
    }
    if(printstats){
      if(!singletonstatsset) singletonstats=new ofstream((outfile+".singleton.stats").c_str(),ofstream::out);
      if(!truncatedstatsset) truncatedstats=new ofstream((outfile+".pair.stats").c_str(),ofstream::out);
    }
  }


  if(collapse) calc_error_matrix();
  
  // Set mismatch threshold
  if(globalmm > 1) globalmm=1.0/globalmm;
  if(usefile2 && globalmm<0)
    globalmm=pairmm;
  if(!usefile2 && globalmm<0)
    globalmm=singlemm;

  // Open file(s), read them line by line, perform alignment of read1 and RC-read2
  // If two files they should have the same number of lines (4 per read pair).

  string read1,read2,temp1,temp2,id1,id2;
  int pairs=0;
  
  istream *readfile1;
  istream *readfile2;
  if(usefile1) readfile1=new ifstream(file1);
  else{
    // Use STDIN
    cerr<<"Expecting data from stdin.\n";
    readfile1=&cin;
  }

  if(usefile2) 
    readfile2=new ifstream(file2);
  else{
    // Use the same read in all alignments constructed from single end adapter
    id2="> Adapter";
    read2=PCR1;
    temp2=string(read2.length(),char(qualitybase+41)); // Create qualities
  }

 // Output some stats and info
  (*settings)<<"Running "<<VERSION<<" using the following options:\n";
 if(usefile2){
   (*settings)<<"Paired end mode\n";
   (*settings)<<"PCR1: "<<PCR1<<((pcr1set)?" (supplied by user)\n":"\n");
   (*settings)<<"PCR2: "<<PCR2<<((pcr2set)?" (supplied by user)\n":"\n");
 }
 else{ 
   (*settings)<<"Single end mode\n";
   (*settings)<<"PCR1: "<<PCR1<<((pcr1set)?" (supplied by user)\n":"\n");
 }
 if(fiveprimeset) (*settings)<<"Trimming 5' end for "<<FIVEPRIME<<"\n";
 (*settings)<<"Alignment shift value: "<<shift<<"\n";
 (*settings)<<"Global mismatch threshold: "<<globalmm<<"\n";
 (*settings)<<"Quality base: "<<qualitybase<<"\n";
 (*settings)<<"Trimming Ns: "<<((removeNs)?"Yes":"No")<<"\n";
 (*settings)<<"Trimming "<<minqual<<": "<<((trimquals)?"Yes":"No")<<"\n";
 (*settings)<<"Minimum genomic length: "<<minimumlength<<"\n";
 (*settings)<<"Collapse overlapping reads: "<<((collapse)?"Yes":"No")<<"\n";
 (*settings)<<"Minimum overlap (in case of collapse): "<<minalignmentlength<<"\n";

  if(DEBUG) cerr<<"TEST\n";
  bool arewegood=false;  
  arewegood=(*readfile1).good();
  if(usefile2) arewegood=(arewegood && (*readfile2).good());

  if(arewegood){
    int counter=0; // We start two lines before the first read
    while (! (*readfile1).eof() && (!usefile2 || !(*readfile2).eof()))
    {
      getline((*readfile1),temp1);
      if(usefile2) getline((*readfile2),temp2);
      counter++;
      if(counter==1){
	id1=temp1;
	if(usefile2) id2=temp2;
      }
      else if(counter==2){
	// Sequences
	read1=toupper(temp1);
	if(usefile2) read2=toupper(temp2);
	pairs++;
      }
      else if(counter==4){
	// Qualities
	lines++;
	counter=0;
	if(fiveprimeset) trim5primeend(read1,temp1);
	pairwisealignment(read1,read2,temp1,temp2,id1,id2,*output1,*output2,*discarded,*singleton,*collapsed,*collapsedtruncated,*singletonstats,*truncatedstats);
      }
    }
  }
  else{
    (*settings) << "Problem opening files.\n\n"; 
    return 1;
  }

  (*settings)<<"\nTotal number of "<<((usefile2)?"read pairs: ":"reads: ")<<lines<<"\n";
  (*settings)<<"Number of unaligned "<<((usefile2)?"read pairs: ":"reads: ")<<noalignment<<"\n";
  (*settings)<<"Number of well aligned "<<((usefile2)?"pairs: ":"reads: ")<<alignedok<<"\n";
  (*settings)<<"Number of inadequate alignments: "<<alignedcrappy<<"\n";
  (*settings)<<"Number of discarded mate 1 reads: "<<discard1<<"\n";
  (*settings)<<"Number of singleton mate 1 reads: "<<keep1<<"\n";
  if(usefile2) (*settings)<<"Number of discarded mate 2 reads: "<<discard2<<"\n";
  if(usefile2) (*settings)<<"Number of singleton mate 2 reads: "<<keep2<<"\n";

  (*settings)<<"\nNumber of "<<((usefile2)?"read pairs":"reads")<<" with adapter: "<<numberofseqswithadapter<<"\n";
  if(collapse){
    (*settings)<<"Number of full-length collapsed pairs: "<<numberoffulllengthcollapsed<<"\n";
    (*settings)<<"Number of truncated collapsed pairs: "<<numberoftruncatedcollapsed<<"\n";
  }
  (*settings)<<"Number of retained reads: "<<totalnumberofgoodreads<<"\n";
  (*settings)<<"Number of retained nucleotides: "<<totalnumberofnucleotides<<"\n";		  
  (*settings)<<"Average read length of trimmed reads: "<<((totalnumberofgoodreads>0)?((double) totalnumberofnucleotides/totalnumberofgoodreads):0)<<"\n";

  // Close files

  // Always open
  (*discarded).flush();
  (*settings).flush();
  (*output1).flush();

  // Always open if stats are printed
  if(printstats){
    (*truncatedstats).flush();
  }

  // Open if paired end
  if(usefile2){
    (*output2).flush();
    (*singleton).flush();
    // Open if PE and collapse
    if(collapse){
      (*collapsed).flush();
      (*collapsedtruncated).flush();
    }
    // Open if PE and stats are printed
    if(printstats){
      (*singletonstats).flush();
    }
  }
  return 0;
}
