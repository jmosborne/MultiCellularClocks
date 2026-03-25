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

#ifndef CIRCADIANREPULSIONFORCE_HPP_
#define CIRCADIANREPULSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "RepulsionForce.hpp"

/**
 * Repulsion force with thresholded circadian modulation.
 *
 * The usual RepulsionForce law is retained, but the pairwise spring constant
 * is multiplied by a threshold function of the mean circadian cycle value of
 * the two interacting cells:
 *
 *   multiplier = 0         if mean_circadian_cycle < threshold
 *   multiplier = amplitude if mean_circadian_cycle >= threshold
 *
 * Requires CircadianRhythmModifier to populate "circadian_cycle" in CellData.
 */
template<unsigned DIM>
class CircadianRepulsionForce : public RepulsionForce<DIM>
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<RepulsionForce<DIM> >(*this);
        archive & mAmplitude;
        archive & mThreshold;
    }

    /** Spring multiplier used when the pair is above threshold. */
    double mAmplitude;

    /** Threshold applied to the mean circadian_cycle of a node pair. */
    double mThreshold;

public:
    CircadianRepulsionForce()
        : RepulsionForce<DIM>(),
          mAmplitude(1.0),
          mThreshold(0.0)
    {
    }

    virtual ~CircadianRepulsionForce()
    {
    }

    double GetAmplitude() const
    {
        return mAmplitude;
    }

    void SetAmplitude(double amplitude)
    {
        mAmplitude = amplitude;
    }

    double GetThreshold() const
    {
        return mThreshold;
    }

    void SetThreshold(double threshold)
    {
        mThreshold = threshold;
    }

    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength) override
    {
        double base_multiplier = RepulsionForce<DIM>::VariableSpringConstantMultiplicationFactor(
            nodeAGlobalIndex,
            nodeBGlobalIndex,
            rCellPopulation,
            isCloserThanRestLength);

        CellPtr p_cell_a = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
        CellPtr p_cell_b = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

        double circadian_cycle_a = 0.0;
        double circadian_cycle_b = 0.0;

        if (p_cell_a->GetCellData()->HasItem("circadian_cycle"))
        {
            circadian_cycle_a = p_cell_a->GetCellData()->GetItem("circadian_cycle");
        }
        if (p_cell_b->GetCellData()->HasItem("circadian_cycle"))
        {
            circadian_cycle_b = p_cell_b->GetCellData()->GetItem("circadian_cycle");
        }

        double mean_circadian_cycle = 0.5 * (circadian_cycle_a + circadian_cycle_b);
        double circadian_multiplier = (mean_circadian_cycle < mThreshold) ? 0.0 : mAmplitude;

        return base_multiplier * circadian_multiplier;
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<Amplitude>" << mAmplitude << "</Amplitude>\n";
        *rParamsFile << "\t\t\t<Threshold>" << mThreshold << "</Threshold>\n";

        RepulsionForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircadianRepulsionForce)

#endif // CIRCADIANREPULSIONFORCE_HPP_
