//============================================================================
// Name        : Supernmotifs program
// Author      : Jean-Pierre Sehi Glouzon
// Copyright   : GNU/GPL
// Description : Supernmotifs algorithm in C++, Ansi-style
//============================================================================

#include<string>
#include<stdio.h>
#include <iostream>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctype.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <queue>
#include <stdexcept>
#include "RepMotif.h"
#include "RepWeightedMotif.h"
#include "Eigen/Dense"
#include <cstdio>

using namespace std;

string i_PathFileInput;//input file
string o_PathDirOutput;//Output directory
int p_outputOption;// Output options
int maxLevelNmotifs;//level of nmotifs
float minOcc;//minimum occurrence
string filenameTemp;

void initParameters(string&, string&, int&,int&,float&);
void setParameters(int, char*[] ,string&, string&,int&,int&,float&);

void writeMatALLMotifWeighted(const vector<string>,const Eigen::ArrayXXf, map<string,float>,const string ,const string );
void writeMatNmotifsPositionsInSS(const vector<string> , vector< map<string,vector<int>> > , map<string,float>  ,const string , const string );
string help();

int main(int argc, char *argv[])
{

    initParameters(i_PathFileInput, o_PathDirOutput,p_outputOption,maxLevelNmotifs,minOcc);
    setParameters(argc, argv, i_PathFileInput, o_PathDirOutput, p_outputOption, maxLevelNmotifs,minOcc);

    cout<<endl<<"Running the super n-motifs program..."<<endl;

    clock_t t1,t2,t3,t4;
    t1=clock();t2=clock();t3=clock();t4=clock();

    //Extract n-motifs.
    cout<<endl<<"Extract n-motifs..."<<endl;
    RepNmotif repnMotif(i_PathFileInput,maxLevelNmotifs);

    t1=clock()-t1;
    cout << "Took:" << ((float)t1)/CLOCKS_PER_SEC << " sec."<<endl;

    //Filter and weight n-motifs : n-motifs representation.
    cout<<"Filter and weight n-motifs to build the n-motifs representation ..."<<endl;
    RepWeightedMotif repFilterWeightNmotif(repnMotif.getNmotifsAllStructure(),repnMotif.getNmotifsForEachStructure(),minOcc);
    t2=clock()-t2;
    cout << "Took:" << ((float)t2)/CLOCKS_PER_SEC << " sec."<<endl;

    t3=clock()-t3;
    cout << "Took:" << ((float)t3)/CLOCKS_PER_SEC << " sec."<<endl;

    cout<<endl<<"Compute and Write output matNmRep_SSbyNm.csv & matnmPos.csv files:"<<endl;

    cout<<"SS*n-motifs matrix that is the n-motif representation of SS (matNmRep_SSbyNm.csv) ..."<<endl;
    filenameTemp="matNmRep_SSbyNm.csv";
    writeMatALLMotifWeighted(repnMotif.getHeaders(), repFilterWeightNmotif.getMatAllMotifWeigthed(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(), o_PathDirOutput,filenameTemp);

    //write the nucleotide position associates with n-motifs in the n-motifs representation of SS
    cout<<"n-motifsPosition matrix (matnmPos.csv) ..."<<endl;
    filenameTemp="matnmPos.csv";
    writeMatNmotifsPositionsInSS(repnMotif.getHeaders(),repnMotif.getNmotifsForEachStructureWithPosNucOfnmotifs(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() ,o_PathDirOutput, filenameTemp);

    cout<<"Execution of n-motifs program completed"<<endl;
     t4=clock()-t4;
    cout << "Took:" << ((float)t4)/CLOCKS_PER_SEC << " sec."<<endl;
    return 0;
}

/******************************************************************************/
void initParameters(string& i_PathFileInput,string& o_PathDirOutput, int& p_outputOption,int& maxLevelNmotifs, float& minOcc){
    i_PathFileInput=""; //input file
    o_PathDirOutput=""; //Output directory
    p_outputOption=0; //Output the dissimilarity matrix
    maxLevelNmotifs=1; //n-motifs parameters
    minOcc=0;//automatic minimum support
}

/******************************************************************************/
void setParameters(int argc,char* argv[],string& i_PathFileInput,string& o_PathDirOutput,int& p_outputOption, int& maxLevelNmotifs, float& minOcc){

    ifstream viennaFile;
    //string filename="helpTerminal.txt";
    //string filename="README.md";

    //ifstream helpFile (filename.c_str());

    string lineHelpFile;
    bool paramRequiredinput=false;
    bool paramRequiredoutput=false;

    struct stat sb;
    for (int i=0;i<argc;i++){
        if (argv[i][0]=='-'){
            switch ( argv[i][1] ) {
            case 'i':
                if (argv[i+1]!=nullptr)
                {
                    i_PathFileInput=string(argv[i+1]);
                    viennaFile.open(i_PathFileInput.c_str());
                      if (!(viennaFile.is_open()))
                      {throw invalid_argument("Unable to open vienna file. Please check the input path file.");}
                    viennaFile.close();
                    paramRequiredinput=true;
                }
                else
                {   throw invalid_argument("Empty value for parameter -i. Please enter the input path file.");}
            break;
            case 'o':
                if (argv[i+1]!=nullptr)
                {
                 o_PathDirOutput=string(argv[i+1]);
                 if (!(stat(o_PathDirOutput.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
                  {
                    cerr<<"The ouput path directory doesn't exits. It will be created."<<endl;
                    #if defined(_WIN32)
                        mkdir(o_PathDirOutput.c_str());
                         #else
                        mkdir(o_PathDirOutput.c_str(), 0700);
                         #endif
                  }
                 paramRequiredoutput=true;
                }
                else
                {   throw invalid_argument("Empty value for parameter -o. Please enter the output path directory.");}
              break;

            case 'n':
                if (argv[i+1]!=nullptr)
                {
                    maxLevelNmotifs=stoi(argv[i+1]);
                    if (!(maxLevelNmotifs==0)&&!(maxLevelNmotifs==1)&& !(maxLevelNmotifs==2))
                    {throw invalid_argument("Maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2. "
                                            "For instance: when it is set to 0, 0-motifs will be extracted. "
                                            "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -n. "
                                           "The maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2."
                                           "For instance: when it is set to 0, 0-motifs will be extracted. "
                                           "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");};
              break;

            case 'm':
                if (argv[i+1]!=nullptr)
                {
                    minOcc=stod(argv[i+1]);
                    if (!(minOcc>=0))
                    {throw invalid_argument("Minimum occurrence of motifs parameter (-m) must be positive (>=0)."
                                            "N-motifs with occurrence below the minimum occurrence are removed."
                                            "When it is set to 0, it computes the automatic n-motif minimum occurrence.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -m."
                                           "Minimum occurrence of motifs parameter (-m) must be positive (>=0). "
                                           "N-motifs with occurrence below the minimum occurrence are removed."
                                           "When it is set to 0, it computes the automatic n-motif minimum occurrence.");};

              break;
            case 'h':
                  cout<< help();

                exit(EXIT_SUCCESS);
              break;
            default:
                cerr<<"No parameters found or wrong parameters. Please check parameter list in the help (-h)."<<endl;
                exit(EXIT_FAILURE);
            break;
            }
            if ((argv[i][1]!='h')&&(argv[i][1]!='i')&&(argv[i][1]!='o')&&(argv[i][1]!='p')&&(argv[i][1]!='n')&&(argv[i][1]!='s')&&(argv[i][1]!='m')&&(argv[i][1]!='g'))
            {cerr<<"Wrong parameters. Please check parameter list in the help (-h).";exit(EXIT_FAILURE);}
        }
    }
    if((paramRequiredinput==false || paramRequiredoutput==false))
        {cerr<<"The program requires the input vienna file (-i) and the ouput folder (-o)."<<endl;exit(EXIT_FAILURE);}
}

/******************************************************************************/
void writeMatALLMotifWeighted(const vector<string> headers, const Eigen::ArrayXXf matALLMotifWeighted, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    map<string,float>::iterator it;
    outputMatMotifs<<" ";
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {outputMatMotifs<<"," << it->first;}
    outputMatMotifs<<"\n";

    //write data
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    for (unsigned int i=0;i<matALLMotifWeighted.rows();i++)
    {outputMatMotifs<<headers[i]<<","<<matALLMotifWeighted.row(i).format(Comma)<<endl;}

    outputMatMotifs.close();
}

/******************************************************************************/
void writeMatNmotifsPositionsInSS(const vector<string> headers, vector< map<string,vector<int>> > FeatureForEachStructureWithPosNuc, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){


    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    map<string,float>::iterator it;
    outputMatMotifs<<" ";
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {outputMatMotifs<<"," << it->first;}
    outputMatMotifs<<"\n";

    //write data
    map<string,vector<int>>::iterator it2;
    for(unsigned int i=0;i<headers.size();i++)
    {
        outputMatMotifs<<headers[i]<<",";
        for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
        {

            it2=FeatureForEachStructureWithPosNuc[i].find(it->first);
            if (it2!=FeatureForEachStructureWithPosNuc[i].end())
            {
                for(unsigned int j=0;j<it2->second.size();j++)
                {
                    outputMatMotifs<<it2->second[j]+1<<"|";
                }
            }
            else
            {
                outputMatMotifs<<"x";
            }
            outputMatMotifs<<",";
        }
        outputMatMotifs<<"\n";
    }

    outputMatMotifs.close();
}

/************************************************************************************/

string help()
{
    //From README.md
    string help;

    help=
"\n### n-motifs model for representing RNA secondary structures ###\n"
"Usage:\n"
"\n"
"* nmotifs [Parameters]...\n"
"\n"
"Examples :\n"
"\n"
"* Execute n-motifs in command line\n"
"\n"
"  ./pathOfSupernmotifsProgram/nmotifs -i /PathtoDbFile\n"
"  -o /OutputDirectoryPath/\n"
"\n"
"Important notes:\n"
"\n"
"* Circular RNA\n"
"\n"
"  For the processing of secondary structures of circular RNA, adding 'c_'\n"
"   at the beginning of the header of each circular RNA is required.\n"
"\n"
"* Pseudoknots and Gquadruplexes\n"
"\n"
"  Special characters '{}','<>','[]', and alphabets such as 'Aa','Bb','Zz'\n"
"  are used to represent base pairs involved in pseudoknots. '+' is used to\n"
"  represent each guanine involved in the Gquadraplexes formation.\n"
"\n"
"### Parameters ###\n"
"\n"
"  -h\n"
"\n"
"	Print help and exit the program.\n"
"\n"
"  -i [input vienna file]\n"
"\n"
"    Input file of RSS in vienna/db format (required).\n"
"    >strucID\n"
"    AAAAAUU\n"
"    ((...))\n"
"\n"
"  -o [output folder]\n"
"\n"
"    Results folder (required).\n"
"\n"
"  -n [0|1|2]\n"
"\n"
"    Specify the maximum level of n-motifs used to\n"
"    extract the n-motifs. It must be 0, 1 or 2. For instance: when\n"
"    it is set to 0, 0-motifs will be extracted. When it is set\n"
"    to 1, 0-motifs and 1-motifs will be extracted. etc.\n"
"    (default : -n 1)\n"
"\n"
"  -m [m>=0]\n"
"\n"
"    Specify the minimum occurrence of n-motifs. N-motifs with occurrence\n"
"    below the minimum occurrence are removed. When it is set to 0,\n"
"    it computes the automatic n-motif minimum occurrence.\n"
"    (default : -m 0)\n"
"\n";

return help;
}
