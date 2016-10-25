#include "solarsystem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#define M_PI       3.14159265358979323846

using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_maxDistance(0.0),
    m_minDistance(1e+25),
    m_generalRelativity(false),
    m_timesteps(0)
{
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass, double radius, string name, int ID) {
    m_bodies.push_back( CelestialBody(position, velocity, mass, radius, name, ID) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy() {
    /*This function calucaltes the forces and energy between every pair of celestial bodies
     *and updates the forces and that are acting up each celestial body.
    */

    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            //Forces are given by: F_1 = - F_2 = G*(M_1*M_2)/(r*r).
            //M_1 is the mass of the first celestial body, and M_2 is of the second.
            //G is the gravitational constant = 6.674×10−11 m^3*kg^{-1}s^{-2}*
            //F_1 is the force of the first celestial body and F_2 is of the second.
            //r is the distance between the centre of masses of the two celestial bodies.
            double G = 4*M_PI*M_PI; double c = 63239.7263; //This is lightspeed in AU/yr.
            int gR;
            if(m_generalRelativity) { //Here the generalRelativity factor is either added or not.
                gR = 1;
            } else {
                gR = 0;
            }
            double l = (deltaRVector.cross(body2.velocity)).length();
            vec3 F_1 = ((G*body1.mass*body2.mass/(dr*dr*dr))*(1.0+gR*((3.0*l*l)/(dr*dr*c*c))))*deltaRVector;
            body1.force -= F_1;
            body2.force += F_1;

            m_potentialEnergy += -G*body1.mass*body2.mass/dr;
        }

        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
        m_angularMomentum += body1.position.cross(body1.mass*body1.velocity);


    }
    m_timesteps++;
}

int SolarSystem::numberOfBodies() const {
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const {
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const {
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const {
    return m_kineticEnergy;
}

bool SolarSystem::collisionDetection(CelestialBody body1, CelestialBody body2) const {
    /*This function detects if two celestial bodies have crashed into eachother.
    */
    bool collision = false;
    vec3 deltaRVector = body1.position - body2.position;
    double dr = deltaRVector.length();
    if(dr <= body1.radius + body2.radius) {
        collision = true;
    }
    return collision;
}

void SolarSystem::writeToFile(string filename) {
    if(!m_ofile.good()) {
        m_ofile.open(filename.c_str(), ofstream::out);
        if(!m_ofile.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_ofile << numberOfBodies() << endl;
    m_ofile << "Comment line that needs to be here. Balle." << endl;
    //This ID is for Ovito program to give the planets different sizes and colours.
    for(CelestialBody &body : m_bodies) {
        m_ofile << body.ID << " " << body.radius << " " << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n";
    }
}

void SolarSystem::writeToFileLogarithm(string filename){
    //In this function we are going to scale the output back from AU to km.
    //And after that we are going to find the polar coordinates and find the logarithm of the radius
    //such that the distances between the celestial bodies are logarithmic.
    //Unfortunatly the result is bizarre.
    if(!m_ofile.good()) {
        m_ofile.open(filename.c_str(), ofstream::out);
        if(!m_ofile.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }
    m_ofile << numberOfBodies() << endl;
    m_ofile << "Comment line that needs to be here. Balle." << endl;
    double conversion = 149597871; //1 AU = 149597871 km
    //x = r*sin(theta)*cos(phi)
    //y = r*sin(theta)*sin(phi)
    //z = r*cos(theta)
    double phi; double theta; double r; double rLog;
    vec3 convertedPosition;
    for(CelestialBody &body : m_bodies) {
        convertedPosition = body.position*conversion;
        if(convertedPosition.x() != 0 || convertedPosition.y() != 0 || convertedPosition.z() != 0) {
            r = convertedPosition.length();
            if(convertedPosition.x() != 0 && convertedPosition.y() != 0) {
                phi = atan2(convertedPosition.y(),convertedPosition.x()) + M_PI;
                theta = acos(convertedPosition.z()/r);
                rLog = log(r);
                convertedPosition.setX(rLog*sin(theta)*cos(phi));
                convertedPosition.setY(rLog*sin(theta)*sin(phi));
                convertedPosition.setZ(rLog*cos(theta));
            }
        }

        m_ofile << body.ID << " " <<  convertedPosition.x() << " " << convertedPosition.y() << " " << convertedPosition.z() << "\n";

        convertedPosition.zeros();
    }
}

void SolarSystem::findMaxMinDistance(string body1Name, string body2Name) {
    /*This function finds the maximum and minimum distance between two given celestialBodies.
     *It will then update m_maxDistance and m_minDistance so that at the end of the run
     *only the minimum and maximum values will be stored.
    */


    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            if(body1.name.compare(body1Name) == 0 && body2.name.compare(body2Name) == 0) {
                if(m_minDistance >= deltaRVector.length()) {
                    m_minDistance = deltaRVector.length();
                    setPerihelionCoordinates(deltaRVector);
                }
                if(m_maxDistance <= deltaRVector.length()) {
                    m_maxDistance = deltaRVector.length();
                }
            } else if (body1.name.compare(body2Name) == 0 && body2.name.compare(body1Name) == 0) {
                if(m_minDistance >= deltaRVector.length()) {
                    m_minDistance = deltaRVector.length();
                    setPerihelionCoordinates(deltaRVector);
                }
                if(m_maxDistance <= deltaRVector.length()) {
                    m_maxDistance = deltaRVector.length();
                }
            }
        }
    }
}

double SolarSystem::maxDistance(){
    return m_maxDistance;
}

double SolarSystem::minDistance() {
    return m_minDistance;
}

void SolarSystem::findCoordinates(CelestialBody body1, CelestialBody body2, double minDist) {
    /*This function attempts to find coordinates of the further away object that is 'minDist' away,
     * and set the m_perihelionCoordinates to the coordinates found here.
     * If Body1 is at (0,0,0) and Body2 is at (1,0,0) and minDist is 1. Then this program set
     *m_periihelionCoordinates to (1,0,0).
    */

    double tolerance = 1e-6;
    vec3 deltaRVector = body1.position - body2.position;
    if(abs(deltaRVector.length() - minDist) <= tolerance) {
        if(body1.position.length() < body2.position.length()) {
            setPerihelionCoordinates(body2.position);
            cout << body2.position << endl;
            setTheta(atan2(m_perihelionCoordinates.y(),m_perihelionCoordinates.x()));
        } else {
            setPerihelionCoordinates(body1.position);
            cout << body1.position << endl;
            setTheta(atan2(m_perihelionCoordinates.y(),m_perihelionCoordinates.x()));
        }
    }
}

void SolarSystem::setTheta(double theta) {
     m_theta.push_back(theta);
}

vector<double> SolarSystem::getTheta() const{
    return m_theta;
}

void SolarSystem::setPerihelionCoordinates(vec3 perihelionCoordinates) {
    m_perihelionCoordinates = perihelionCoordinates;
}

vec3 SolarSystem::getPerihelionCoordinates() const{
    return m_perihelionCoordinates;
}

void SolarSystem::findPerihelion(CelestialBody body1, CelestialBody body2) {
    if(m_timesteps == 2) {
        m_prevPrevCoord = body2.position-body1.position;
    }
    if(m_timesteps == 3) {
        m_prevCoord = body2.position-body1.position;
    }
    if(m_timesteps > 3) {
        m_currentCoord = body2.position-body1.position;

        if((m_prevPrevCoord.length() > m_prevCoord.length()) && (m_currentCoord.length() > m_prevCoord.length())) {
            m_periCoord.push_back(m_prevCoord);
            setTheta(atan2(m_periCoord.back().y(),m_periCoord.back().x()));
        }
        m_prevPrevCoord = m_prevCoord;
        m_prevCoord = m_currentCoord;

    }
}

vector<vec3> SolarSystem::getPeriCoord() {
    return m_periCoord;
}

void SolarSystem::setGeneralRelativity() {
    m_generalRelativity = true;
}

vec3 SolarSystem::angularMomentum() const {
    return m_angularMomentum;
}

void SolarSystem::setCenterOfMass() {
    findCenterOfMass();
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body = m_bodies[i];
        body.position -= m_centerOfMass;
    }
}

void SolarSystem::findCenterOfMass() {
    /*This function attempts to find the center of mass of the solarsystem
     * and set m_centerOfMass to this value.
    */
    m_centerOfMass.zeros();
    double sumOfMass = 0;
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body = m_bodies[i];
        sumOfMass += body.mass;
        m_centerOfMass += body.mass*body.position;
    }
    m_centerOfMass = m_centerOfMass/sumOfMass;
}

vec3 SolarSystem::getCenterOfMass() const{
    return m_centerOfMass;
}

void SolarSystem::setMomentum() {
    /*This function attempts to set the velocity of the first celestialBody (Usually the sun)
     *so that the mass momentum is conserved and the center of mass does not wobble.
    */
    m_momentum.zeros();
    for(int i=1; i<numberOfBodies(); i++) {
        CelestialBody &body = m_bodies[i];
        m_momentum += body.mass*body.velocity;
    }
    CelestialBody &centralCelestialObject = m_bodies[0];
    centralCelestialObject.velocity = (-1)*m_momentum/centralCelestialObject.mass;
}

vec3 SolarSystem::getMomentum() const{
    return m_momentum;
}

std::vector<CelestialBody> &SolarSystem::bodies() {
    return m_bodies;
}

ofstream &SolarSystem::file()
{
    return m_ofile;
}
