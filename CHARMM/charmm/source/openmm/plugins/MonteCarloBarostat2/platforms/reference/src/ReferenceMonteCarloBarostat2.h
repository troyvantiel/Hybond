/* Portions copyright (c) 2010 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#ifndef __ReferenceMonteCarloBarostat2_H__
#define __ReferenceMonteCarloBarostat2_H__

#include <openmm/reference/RealVec.h>
#include <vector>

class ReferenceMonteCarloBarostat2 {
private:
    std::vector<RealOpenMM> savedAtomPositions[3];
    std::vector<std::vector<int> > molecules;

public:
    ReferenceMonteCarloBarostat2(int numAtoms, const std::vector<std::vector<int> >& molecules);

    ~ReferenceMonteCarloBarostat2();

    /**
     * Apply the barostat at the start of a time step.
     *
     * @param atomPositions      atom positions
     * @param boxSize            the periodic box dimensions
     * @param scale              the factor by which to scale atom positions
     */
    void applyBarostat(std::vector<OpenMM::RealVec>& atomPositions, const OpenMM::RealVec& boxSize, OpenMM::Vec3& scale);

    /**
     * Restore atom positions to what they were before applyBarostat() was called.
     *
     * @param atomPositions      atom positions
     */
    void restorePositions(std::vector<OpenMM::RealVec>& atomPositions);
};

#endif
