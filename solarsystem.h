#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>

class SolarSystem
{
public:
    SolarSystem();
    CelestialBody &createCelestialBody(vec3 position, vec3 velocity, double mass, double radius, string name, int ID);
    void calculateForcesAndEnergy();
    int numberOfBodies() const; //const means that whoever calls the function 'numberOfBodies()' will not change the value.
                                //It is good for get-functions where you just want the value, and not to change it, compared
                                //to set-functions.
    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    bool collisionDetection(CelestialBody body1, CelestialBody body2) const;
    void writeToFile(std::string filename);
    void writeToFileLogarithm(std::string filename);
    void findMaxMinDistance(string body1Name, string body2Name);
    double maxDistance();
    double minDistance();
    void setGeneralRelativity();
    vec3 angularMomentum() const;
    void setCenterOfMass();
    void findCenterOfMass();
    vec3 getCenterOfMass() const;
    void setMomentum();
    vec3 getMomentum() const;
    void setPerihelionCoordinates(vec3 perihelionCoordinates);
    vec3 getPerihelionCoordinates() const;
    void findCoordinates(CelestialBody body1Name, CelestialBody body2Name, double minDistance);
    void setTheta(double theta);
    void findPerihelion(CelestialBody body1, CelestialBody body2);
    vector<vec3> getPeriCoord();
    vector<double> getTheta() const;
    void findPerihelion(string body1Name, string body2Name) const;
    std::vector<CelestialBody> &bodies();
    ofstream &file();

private:
    std::vector<CelestialBody> m_bodies; //All the planets, astroids etc.
    vec3 m_angularMomentum; //The solarsystems angular momentum
    std::ofstream m_ofile; //m_ stands for member, member variable, a member of the class. It is just a convention.
    vec3 m_centerOfMass;
    vec3 m_momentum;
    vec3 m_perihelionCoordinates;
    double m_maxDistance;
    double m_minDistance;
    bool m_generalRelativity;
    double m_kineticEnergy;
    double m_potentialEnergy;
    vector<double> m_theta;
    vec3 m_prevPrevCoord;
    vec3 m_prevCoord;
    vec3 m_currentCoord;
    int m_timesteps;
    vector<vec3> m_periCoord;
};

#endif // SOLARSYSTEM_H
