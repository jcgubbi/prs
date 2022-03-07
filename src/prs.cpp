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

const int num_files = 2;

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

map<string, int> getReferences(string a) {
    vector<string> split;
    map<string, int> m;

    while(a.find('\t') != string::npos) {
        size_t pos = a.find('\t');
        split.push_back(a.substr(0,pos));
        a = a.substr(pos+1, a.length() - pos);
    }
    split.push_back(a);

    for(int i = 0; i < split.size(); i++){
        if(split[i] == "rsID")
            m.insert(pair<string, int>("rsID", i));
        if(split[i] == "effect_allele")
            m.insert(pair<string, int>("effect_allele", i));
        if(split[i] == "other_allele")
            m.insert(pair<string, int>("other_allele", i));
        if(split[i] == "effect_weight")
            m.insert(pair<string, int>("effect_weight", i));
        if(split[i] == "allelefrequency_effect")
            m.insert(pair<string, int>("allelefrequency_effect", i));
        if(split[i] == "locus_name")
            m.insert(pair<string, int>("locus_name", i));
        if(split[i] == "OR")
            m.insert(pair<string, int>("OR", i));
    }


    for(map<string, int>::const_iterator it = m.begin();
        it != m.end(); ++it)
    {
        std::cout << it->first << " " << it->second << "\n";
    }

    return m;
}
Line process_line(string a, map<string, int> ref) {
    Line l;

    vector<string> split;

    // we know format of the line
    // tab separated with 7 fields
    while(a.find('\t') != string::npos) {
        size_t pos = a.find('\t');
        split.push_back(a.substr(0,pos));
        a = a.substr(pos+1, a.length() - pos);
    }
    split.push_back(a);

    int rsid_i = ref.find("rsID") != ref.end() ? ref.find("rsID")->second : -1;
    int ef_i = ref.find("effect_allele") != ref.end() ? ref.find("effect_allele")->second : -1;
    int ot_i = ref.find("other_allele") != ref.end() ? ref.find("other_allele")->second : -1;
    int ew_i = ref.find("effect_weight") != ref.end() ? ref.find("effect_weight")->second : -1;
    int af_i = ref.find("allelefrequency_effect") != ref.end() ? ref.find("allelefrequency_effect")->second : -1;
    int ln_i = ref.find("locus_name") != ref.end() ? ref.find("locus_name")->second : -1;
    int or_i = ref.find("OR") != ref.end() ? ref.find("OR")->second : -1;

    l.rsid = rsid_i != -1 ? split[rsid_i] : "unknown";
    l.effect = ef_i != -1 ? split[ef_i].at(0) : 'K';
    l.other = ot_i != -1 ? split[ot_i].at(0) : 'K';
    l.weight = ew_i != -1 ? stod(split[ew_i]) : 0;
    l.frequency = af_i != -1 ? stod(split[af_i]) : 0.5;
    l.locus = ln_i != -1 ? split[ln_i] : "unknown";
    l.or_v = or_i != -1 ? stod(split[or_i]): 0;

    return l;
}

Individual randomPerson(vector<Line> locs){
    Individual r;
    
    r.prs = 0;
    //random allele assignment
    for(int i = 0; i < locs.size(); i++){
        int random = rand()%100;
        auto itr = r.alleles.find(locs[i].rsid);
        if(itr != r.alleles.end() && itr->second == locs[i].effect) {
            r.prs += locs[i].weight;
        }
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

int main (int argc, char *argv[]) {  
    string prs_files[num_files] = {"resources/PGS000011.txt", "resources/PGS000010.txt"};
    vector<Line> locs;

    for (int i = 0; i < num_files; i++){
        string line;
        ifstream myfile (prs_files[i]);
        map<string, int> ref;
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                if(line[0] == '#') {
                    continue;
                }
                if(line.find("rsID") != string::npos) {
                    cout << "----Parsing PGS File: " << prs_files[i] << "----" << endl;
                    ref = getReferences(line);
                    continue;
                }
                locs.push_back(process_line(line, ref));
            }
            myfile.close();
        }
    }

    string n_s = argv[1];
    int n = stoi(n_s);
    vector<Individual> group;

    cout << endl << "Number of individuals to rank: " << n_s << endl << endl;

    cout << "----Creating " << n << " random individuals----" << endl << endl;
    for(int i = 0; i < n; i++){
        group.push_back(randomPerson(locs));
    }

    cout << "----Sorting individuals----" << endl << endl; 
    sort(group.begin(), group.end(), comp);

    cout << "----Results----" << endl;
    cout << "Name\tScore" << endl;

    for(int i = 0; i < n; i++){
        cout << group[i].name << '\t' << group[i].prs << endl;
    }

    double theoretical_max = 0;
    for(int i = 0; i < locs.size(); i++) {
        theoretical_max += locs[i].weight;
    }
    cout << "Theoretical maximum: " << theoretical_max << endl;

    return 0;
}