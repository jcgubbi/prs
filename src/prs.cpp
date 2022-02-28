// basic file operations
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <map>

using namespace std;

//general steps
//1. Read files
//2. Store info in some sort of data structure
//3. Randomize input distribution
//4. Compute PRS of inputs
//5. Sort by PRS

const int num_files = 1;

struct Line {
    string rsid;
    char effect;
    char other;
    double weight;
    double frequency;
    string locus;
    double or_v;
};

struct Individual {
    string name;
    map<string, char> alleles;
    double prs;
};

Line process_line(string a) {
    Line l;

    vector<string> split;

    // we know format of the line
    // tab separated with 7 fields
    for (int i = 0; i < 6; i++) {
        size_t pos = a.find('\t');
        split.push_back(a.substr(0,pos));
        a = a.substr(pos+1, a.length() - pos);
    }
    split.push_back(a);

    l.rsid = split[0];
    l.effect = split[1].at(0);
    l.other = split[2].at(0);
    l.weight = stod(split[3]);
    l.frequency = stod(split[4]);
    l.locus = split[5];
    l.or_v = stod(split[6]);

    return l;
}

Individual randomPerson(vector<Line> locs){
    Individual r;
    
    r.prs = 0;
    //random allele assignment
    for(int i = 0; i < locs.size(); i++){
        int random = rand()%100;
        if ((random * 1.0) / 100 > locs[i].frequency) {
            r.alleles.insert(pair<string, char>(locs[i].rsid, locs[i].effect));
            r.prs += locs[i].weight;
        } else {
            r.alleles.insert(pair<string, char>(locs[i].rsid, locs[i].other));
        }
    }

    //random name
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    string tmp_s;
    tmp_s.reserve(20);

    for (int i = 0; i < 20; ++i) {
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    r.name = tmp_s;

    return r;
}

int comp(Individual a, Individual b) {
    return a.prs > b.prs;
}

int main () {  
    string prs_files[num_files] = {"resources/PGS000011.txt"};
    vector<Line> locs;
    for (int i = 0; i < num_files; i++){
        string line;
        ifstream myfile (prs_files[i]);
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                if(line[0] == '#' || line.find("rsID") != string::npos){
                    continue;
                }
                locs.push_back(process_line(line));
            }
            myfile.close();
        }
    }
    cout << "Number of individuals to rank: ";

    string n_s;
    cin >> n_s;
    int n = stoi(n_s);
    vector<Individual> group;

    for(int i = 0; i < n; i++){
        group.push_back(randomPerson(locs));
    }

    cout << group[0].name << " " << group[0].prs << " " << group[0].alleles.at("rs646776") << endl;
    
    sort(group.begin(), group.end(), comp);

    for(int i = 0; i < n; i++){
        cout << group[i].prs << endl;
    }
    return 0;
}