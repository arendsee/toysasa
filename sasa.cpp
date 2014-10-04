// sasa - sovent accssible surface area calculator
// Zebulun Arendsee (arendsee@iastate.edu)
// BCB569 HW3
// October 2014

#include <math.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
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

void print_atom(struct Atom a)
{
    printf("(%f, %f, %f) %s %s %s (%d, %d) %f\n",
           a.pnt.x, a.pnt.y, a.pnt.z,
           a.residue.c_str(),
           a.element.c_str(),
           a.atom_name.c_str(),
           a.serial_id,
           a.aa_id,
           a.radius
          );
}

void print_point(struct XYZ p)
{
    printf("%f %f %f\n", p.x, p.y, p.z);
}

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

vector<struct Atom> load_pdb_file()
{
    map<string, float> radii;
    radii["C"] = 1.7;
    radii["H"] = 1.2;
    radii["O"] = 1.5;
    radii["N"] = 1.6;
    radii["S"] = 1.8;

    vector<Atom> atoms;
    string line;
    while(cin){
        getline(cin, line);
        if(line.substr(0,4) != "ATOM")
            continue;
        struct Atom atom;
        atom.serial_id = atoi(line.substr(6,4).c_str());
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

bool is_too_close(struct Probe p, vector<struct Atom> a)
{
    for(int i = 0; i < a.size(); i++){
        if(dist(p.pnt, a[i].pnt) < p.radius + a[i].radius){
            return true;
        }
    }
    return false; 
}

vector<struct Atom> get_nearby_atoms(Atom a, vector<struct Atom> av, float d)
{
    vector<struct Atom> nearby;
    for(int i = 0; i < av.size(); i++){
        if(dist(a.pnt, av[i].pnt) < d){
            nearby.push_back(av[i]);
        }
    }
    return nearby;
}

struct XYZ get_random_surface_coordinate(struct XYZ p, double radius)
{
    double r1 = ((double) rand() / RAND_MAX);
    double r2 = ((double) rand() / RAND_MAX);

    double xz = 2 * PI * r1;
    double xy = acos(2 * r2 - 1);

    struct XYZ p2;
    p2.x = p.x + radius * cos(xz) * cos(xy);
    p2.y = p.y + radius * sin(xy);
    p2.z = p.z + radius * sin(xz);

    return p2;
}

int main(int argc, char* argv[])
{
    struct Probe p;
    p.radius = 2;
    vector<struct Atom> a = load_pdb_file(); 
    vector<struct Atom> nearby;
    double nearby_radius, solvent_radius;
    int npos;
    int k = 500;
    for(int i = 0; i < a.size(); i++){
        nearby_radius = a[i].radius + 2*p.radius + 2;
        solvent_radius = a[i].radius + p.radius;
        nearby = get_nearby_atoms(a[i], a, nearby_radius); 
        npos = 0;
        for(int j = 0; j < k; j++){
            p.pnt = get_random_surface_coordinate(a[i].pnt, solvent_radius);
            if(!is_too_close(p, nearby)){
                npos++;
            }
        }
        a[i].sasa = float(npos) / k;
        printf("%f %d %d\n", a[i].sasa, npos, k);
    }
}
