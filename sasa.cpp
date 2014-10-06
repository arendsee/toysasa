// sasa - sovent accssible surface area calculator
// Zebulun Arendsee (arendsee@iastate.edu)
// BCB569 HW3
// October 2014
//
// USAGE: sasa.out < myfile.pdb > output.txt

#include <math.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
using namespace std;

const double PI = acos(-1);

struct XYZ { double x,y,z; };


struct Atom
{
    XYZ pnt;
    string residue;
    string element;
    string atom_name;
    int serial_id;
    int aa_id;
    double radius;
    double sasa;
};


struct Probe
{
    XYZ pnt;
    double radius;
};


void trim(string& s)
{
    string::size_type pos = s.find_last_not_of(' ');
    if(pos != string::npos) {
        if (s.length() != pos + 1)
            s.erase(pos + 1);
            pos = s.find_first_not_of(' ');
        if(pos != 0)
            s.erase(0, pos);
    }
    else s="";
}


/*
 * Extract one Atom struct for each atom in the PDB file
 */
vector<struct Atom> load_pdb_file()
{
    map<string, float> radii;
    radii["C"] = 1.7;
    radii["H"] = 1.2;
    radii["O"] = 1.52;
    radii["N"] = 1.55;
    radii["S"] = 1.8;

    vector<Atom> atoms;
    string line;
    while(cin){
        getline(cin, line);
        if(line.substr(0,4) != "ATOM")
            continue;
        struct Atom atom;
        atom.serial_id = atoi(line.substr(6,5).c_str());
        atom.atom_name = line.substr(12,3);
        atom.residue   = line.substr(17,3);
        atom.aa_id     = atoi(line.substr(22,4).c_str());
        atom.pnt.x     = atof(line.substr(30, 8).c_str());
        atom.pnt.y     = atof(line.substr(38, 8).c_str());
        atom.pnt.z     = atof(line.substr(46, 8).c_str());
        atom.element   = line.substr(76,2);
        trim(atom.atom_name);
        trim(atom.residue);
        trim(atom.element);
        atom.radius = radii[atom.element];
        atoms.push_back(atom);
    }
    return atoms;
}


/*
 * Calculate distance between two points
 */
double dist(struct XYZ a, struct XYZ b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + 
                (a.y - b.y) * (a.y - b.y) +
                (a.z - b.z) * (a.z - b.z));
}


/*
 * Deterimines whether a point on the solvenation sphere is aaccessible
 *   TRUE if p does not overlap with any atom in a
 *   FALSE if it does overlap
 */
bool accessible(struct Probe p, vector<struct Atom> a)
{
    for(int i = 0; i < a.size(); i++){
        if(dist(p.pnt, a[i].pnt) < (p.radius + a[i].radius)){
            return false;
        }
    }
    return true; 
}


/*
 * Find all atoms that are within radius d of a given atom
 */
vector<struct Atom> get_nearby_atoms(Atom a, vector<struct Atom> av, float d)
{
    vector<struct Atom> nearby;
    for(int i = 0; i < av.size(); i++){
        if(dist(a.pnt, av[i].pnt) < d){
            // Don't include self in the list of nearby atoms
            if(a.serial_id == av[i].serial_id) {
                continue;
            }
            nearby.push_back(av[i]);
        }
    }
    return nearby;
}


/*
 * Get a random point on a sphere surface centered at p with radius r
 */
struct XYZ get_random_surface_coordinate(struct XYZ p, double r)
{
    double phi = 4 * PI * ((double) rand() / RAND_MAX) - 2 * PI;
    double z = 2 * r * ((double) rand() / RAND_MAX) - r;
    double x = sqrt(r*r - z*z) * cos(phi);
    double y = sqrt(r*r - z*z) * sin(phi);
    struct XYZ p2;
    p2.x = x + p.x;
    p2.y = y + p.y;
    p2.z = z + p.z;
    return p2;
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    // Initialize total van der Waals surface variable
    double total_waals = 0;
    // Initialize solvent accessible surface area variable
    double total_sasa  = 0;
    // Initalize total atomic surface area variable
    double total_area  = 0;

    vector<struct Atom> a = load_pdb_file(); 

    // Initialize probe with default radius
    struct Probe p;
    p.radius = 1.4;

    // Set number of random positions to probe for each atom
    int k = 2000;

    // For each atom, calculate SASA
    for(int i = 0; i < a.size(); i++){

        // Get nearby atoms (this is a speed optimization to avoid calculating
        // distances to every atom in the protein in the monte carlo loop)
        double nearby_radius = a[i].radius + 2*p.radius + 2;
        vector<struct Atom> nearby = get_nearby_atoms(a[i], a, nearby_radius); 

        double solvent_radius = a[i].radius + p.radius;
        int solvent_accessible = 0;

        // Count number of accessible positions on solvent surface
        for(int j = 0; j < k; j++){
            p.pnt = get_random_surface_coordinate(a[i].pnt, solvent_radius);
            if(accessible(p, nearby)){
                solvent_accessible++;
                // printf("%8.3f%8.3f%8.3f  %d  %d\n", p.pnt.x, p.pnt.y, p.pnt.z);
            }
        }


        // Calculate approximate SASA (accessible / total)
        double ratio = float(solvent_accessible) / k;

        double atomic_area = 4 * PI * a[i].radius * a[i].radius;
        double shell_area = 4 * PI * solvent_radius * solvent_radius;

        double sasa = ratio * shell_area;
        double waals = ratio * atomic_area;

        total_waals += waals;
        total_sasa += sasa;
        total_area += atomic_area;
        a[i].sasa = sasa;
    }

    // Print atomic statistics
    printf("Results for each atom\n");
    for(int i = 0; i < a.size(); i++){
        printf("%s-%-2d %3s %8.3f\n", a[i].residue.c_str(), a[i].aa_id, a[i].atom_name.c_str(), a[i].sasa);
    }

    // Print amino acid statistics
    printf("\n");
    printf("Results for each residue\n");
    double total = 0;
    int aa_id = 0;
    string residue;
    for(int i = 0; i < a.size(); i++){
        if(a[i].aa_id != aa_id && aa_id > 0){
            printf("%s-%-2d %8.3f\n", residue.c_str(), aa_id, total);
            total = 0;
        }
        aa_id = a[i].aa_id;
        residue = a[i].residue;
        total += a[i].sasa;
        if(i == a.size() - 1){
            printf("%s-%-2d %8.3f\n", residue.c_str(), aa_id, total);
        }
    }

    // Print totals
    printf("\n");
    printf("Total atomic area: %.1f\n", total_area); 
    printf("Total van der Waals surface: %.1f\n", total_waals); 
    printf("Total solvent accessible surface area: %.f\n", total_sasa); 
}
