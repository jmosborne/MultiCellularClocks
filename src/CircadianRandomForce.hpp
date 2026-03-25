/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CIRCADIANRANDOMFORCE_HPP_
#define CIRCADIANRANDOMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A force that applies random motion to cells, with magnitude modulated by
 * the circadian cycle value.
 *
 * For each cell, a random direction is chosen and a force magnitude is
 * computed as:
 *
 *   F(t) = step_modulation
 *
 * where
 *   step_modulation = 0         if circadian_cycle < threshold
 *   step_modulation = amplitude if circadian_cycle >= threshold
 *
 * where amplitude sets the force magnitude in the active part of the
 * circadian cycle.
 *
 * Requires CircadianRhythmModifier to be running so that
 * "circadian_cycle" is available in cell data.
 */
template<unsigned DIM>
class CircadianRandomForce : public AbstractForce<DIM, DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM, DIM> >(*this);
        archive & mAmplitude;
        archive & mThreshold;
    }

    /** Circadian modulation amplitude (default: 0.5). */
    double mAmplitude;

    /** Circadian threshold for switching modulation (default: 0.0). */
    double mThreshold;

public:

    /** Constructor. */
    CircadianRandomForce();

    /** Destructor. */
    virtual ~CircadianRandomForce();

    /**
     * @return circadian modulation amplitude.
     */
    double GetAmplitude() const;

    /**
     * Set circadian modulation amplitude.
     *
     * @param amplitude modulation amplitude (typically in [0, 1])
     */
    void SetAmplitude(double amplitude);

    /**
     * @return circadian threshold.
     */
    double GetThreshold() const;

    /**
     * Set circadian threshold.
     *
     * @param threshold threshold applied to circadian_cycle
     */
    void SetThreshold(double threshold);

    /**
     * Overridden AddForceContribution() method.
     *
     * Applies a random force to each cell, with magnitude modulated
      * by circadian_cycle from cell data.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircadianRandomForce)

#endif // CIRCADIANRANDOMFORCE_HPP_
