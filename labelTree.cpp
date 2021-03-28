#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;

const int n_pos = 4; // number of positions to be annotated
const string pos[n_pos] = {"222", "439", "477", "614"}; // positions to be annotated in the same order in the csv file
string annt_fileName = "combined_aaAnnotation.csv";
string tree_fileName = "large_with_Australia_tree.txt";

bool replace(std::string& str, const std::string& from, const std::string& to) {
  // obtained from https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

int main()
{
    // Get annotation data from input file
    ifstream annt(annt_fileName.c_str());
    if (!annt.is_open())
    {
        exit(EXIT_FAILURE);
    }
    vector<vector<string>> sequences;
    vector<string> sequence;
    string str, data;
    getline(annt, str); // skip the first line

    while (getline(annt, str))
    {
      str.erase(remove(str.begin(), str.end(), '"'), str.end()); //remove " "
      istringstream iss(str);
      vector<string> sequence;
      for (int i = 0; i < n_pos+1; i++) {
        getline(iss, data, ',');
        sequence.push_back(data);
      }
      sequences.push_back(sequence);
    }
    annt.close();

    //Open the tree file
    ifstream tree_file(tree_fileName.c_str());
    if (!tree_file.is_open())
    {
        exit(EXIT_FAILURE);
    }
    string tree;
    tree_file >> tree;
    tree_file.close();

    // Label of sequences in the tree
    string id, id_label;
    for (int i = 0; i < sequences.size(); i++) {
      id = sequences[i][0];
      id_label = id;
      for (int j = 0; j < n_pos; j++) {
        id_label += ("|" + sequences[i][j+1] + pos[j]);
      }
      replace(tree, id, id_label);
    }
    // Write the labelled tree to a new file
    ofstream fout ("labelled_tree.txt");
    fout << tree;
    fout.close();
    return 0;
}
