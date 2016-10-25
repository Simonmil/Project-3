#include "celestialbody.h"

CelestialBody::CelestialBody(vec3 pos, vec3 vel, double mass_, double radius_, string name_, int ID_) {
    position = pos;
    velocity = vel;
    mass = mass_;
    radius = radius_;
    name = name_;
    ID = ID_;
}

CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass_, double radius_, string name_, int ID_) {
    position = vec3(x,y,z);
    velocity = vec3(vx,vy,vz);
    mass = mass_;
    name = name_;
    radius = radius_;
    ID = ID_;
}

void CelestialBody::resetForce() {
    force.zeros();
}

