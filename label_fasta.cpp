#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
string fastaFile = "gisaid_Denmark_20201109_to_20210208_aligned.fasta";
string pangolinFile = "pangolin_Denmark_20201109_to_20210208.csv";

int main() {
  ifstream fasta(fastaFile.c_str());
  if (fasta.fail()) {
    cout << ".fasta file not found" << endl;
    exit(1);
  }

  ifstream pangolin(pangolinFile.c_str());
  if (pangolin.fail()) {
    cout << "pangolin file not found" << endl;
    exit(1);
  }

  ofstream labelled(("labelled_"+fastaFile).c_str());
  if (labelled.fail()){
    cout << "failed to create the output file" << endl;
  }

  string f_line, p_line, lineage;
  istringstream iss;

  getline(pangolin, p_line); // skip the column names

  while(getline(fasta, f_line)) {
    if (f_line.find(">") == 0) {
      getline(pangolin, p_line);
      iss.str(p_line);
      getline(iss, lineage, ','); // the first variable is not lineage (ID)
      getline(iss, lineage, ',');
      f_line += ("|" + lineage);
    }
    labelled << f_line << endl;
  }

  fasta.close();
  pangolin.close();
  labelled.close();
  return 0;
}
