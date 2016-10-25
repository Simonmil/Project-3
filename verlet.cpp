#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet(double dt) :
    m_dt(dt),
    m_firstStep(true)
{

}

void Verlet::integrateOneStep(SolarSystem &system)
{
    if(m_firstStep) {
        system.calculateForcesAndEnergy();
        m_firstStep = false;
    }

    for(int i=0; i<system.numberOfBodies(); i++) {
        CelestialBody &body1 = system.bodies()[i];
        for(int j=i+1; j<system.numberOfBodies(); j++) {
            CelestialBody &body2 = system.bodies()[j];

            if(body1.mass != 0.0 || body2.mass != 0.0) {
            //If there is a collision. We update the velocity.
                if(system.collisionDetection(body1, body2)) {
                    if(body1.radius > body2.radius) {
                        body1.velocity += body2.velocity*body2.mass/body1.mass;
                        body2.resetForce();
                        body2.radius = 0.0;
                        body2.velocity.zeros();
                        body2.mass = 0.0;
                    } else {
                        body2.velocity += body1.velocity*body1.mass/body2.mass;
                        body1.resetForce();
                        body1.radius = 0.0;
                        body1.velocity.zeros();
                        body1.mass = 0.0;
                    }
                }
            }
        }
    }

    for(CelestialBody &body : system.bodies()) {
        if(body.mass != 0.0){
            body.velocity += 0.5*body.force/body.mass*m_dt;
            body.position += body.velocity*m_dt;
        }
    }

    system.calculateForcesAndEnergy();

    for(CelestialBody &body : system.bodies()) {
        if(body.mass != 0.0){
            body.velocity += 0.5*body.force/body.mass*m_dt;
        }
    }
}
