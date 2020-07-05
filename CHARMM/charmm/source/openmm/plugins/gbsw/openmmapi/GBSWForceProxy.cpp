/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "GBSWForceProxy.h"
#include "GBSWForce.h"

#include "openmm/Force.h"
#include "openmm/serialization/SerializationNode.h"

#include <sstream>

using namespace OpenMM;
using namespace OpenMMGBSW;
using namespace std;

GBSWForceProxy::GBSWForceProxy() : SerializationProxy("GBSWForce") {
}

void GBSWForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const GBSWForce& force = *reinterpret_cast<const GBSWForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setDoubleProperty("soluteDielectric", force.getSoluteDielectric());
    node.setDoubleProperty("solventDielectric", force.getSolventDielectric());
    node.setDoubleProperty("surfaceAreaEnergy", force.getSurfaceAreaEnergy());
    node.setDoubleProperty("rfDielectric", force.getReactionFieldDielectric());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius;
        force.getParticleParameters(i, charge, radius);
        particles.createChildNode("Particle").setDoubleProperty("q", charge).setDoubleProperty("r", radius);
    }
}

void* GBSWForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    GBSWForce* force = new GBSWForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((GBSWForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setSoluteDielectric(node.getDoubleProperty("soluteDielectric"));
        force->setSolventDielectric(node.getDoubleProperty("solventDielectric"));
        force->setReactionFieldDielectric(node.getDoubleProperty("rfDielectric"));
        if (version > 1)
            force->setSurfaceAreaEnergy(node.getDoubleProperty("surfaceAreaEnergy"));
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getDoubleProperty("q"), particle.getDoubleProperty("r"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
