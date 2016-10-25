#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include "vec3.h"
#include <string>

using namespace std;

class CelestialBody
{
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;
    double radius;
    string name;
    int ID;

    CelestialBody(vec3 position, vec3 velocity, double mass, double radius, string name, int ID);
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass, double radius, string name, int ID);
    void resetForce();
};

#endif // CELESTIALBODY_H
