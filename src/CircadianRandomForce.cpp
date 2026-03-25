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

#include "CircadianRandomForce.hpp"

#include <cassert>
#include <cmath>

#include "NonCyclingCellProperty.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
CircadianRandomForce<DIM>::CircadianRandomForce()
    : AbstractForce<DIM, DIM>(),
            mAmplitude(0.5),
            mThreshold(0.0)
{
}

template<unsigned DIM>
CircadianRandomForce<DIM>::~CircadianRandomForce()
{
}

template<unsigned DIM>
double CircadianRandomForce<DIM>::GetAmplitude() const
{
    return mAmplitude;
}

template<unsigned DIM>
void CircadianRandomForce<DIM>::SetAmplitude(double amplitude)
{
    mAmplitude = amplitude;
}

template<unsigned DIM>
double CircadianRandomForce<DIM>::GetThreshold() const
{
    return mThreshold;
}

template<unsigned DIM>
void CircadianRandomForce<DIM>::SetThreshold(double threshold)
{
    mThreshold = threshold;
}

template<unsigned DIM>
void CircadianRandomForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    for (auto p_cell = rCellPopulation.Begin(); p_cell != rCellPopulation.End(); ++p_cell)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*p_cell);
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

        // Get circadian cycle value from cell data, default to 0 if not present
        double circadian_cycle = 0.0;
        if (p_cell->GetCellData()->HasItem("circadian_cycle"))
        {
            circadian_cycle = p_cell->GetCellData()->GetItem("circadian_cycle");
        }

        // Compute force magnitude with thresholded circadian modulation.
        // Non-cycling cells always receive zero random force.
        // Otherwise, magnitude = 0 if circadian_cycle < threshold, else amplitude.
        const double magnitude = p_cell->template HasCellProperty<NonCyclingCellProperty>()
            ? 0.0
            : ((circadian_cycle < mThreshold) ? 0.0 : mAmplitude);

        // Generate a random unit direction using angles.
        c_vector<double, DIM> random_direction = zero_vector<double>(DIM);

        if (DIM == 1)
        {
            random_direction[0] = (RandomNumberGenerator::Instance()->ranf() < 0.5) ? -1.0 : 1.0;
        }
        else if (DIM == 2)
        {
            const double angle = 2.0 * 3.14159265358979323846 * RandomNumberGenerator::Instance()->ranf();
            random_direction[0] = std::cos(angle);
            random_direction[1] = std::sin(angle);
        }
        else if (DIM == 3)
        {
            const double azimuth = 2.0 * 3.14159265358979323846 * RandomNumberGenerator::Instance()->ranf();
            const double cos_polar = 2.0 * RandomNumberGenerator::Instance()->ranf() - 1.0;
            const double sin_polar = std::sqrt(1.0 - cos_polar*cos_polar);

            random_direction[0] = sin_polar * std::cos(azimuth);
            random_direction[1] = sin_polar * std::sin(azimuth);
            random_direction[2] = cos_polar;
        }

        // Apply force
        c_vector<double, DIM> force = magnitude * random_direction;
        p_node->AddAppliedForceContribution(force);
    }
}

template<unsigned DIM>
void CircadianRandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Amplitude>" << mAmplitude << "</Amplitude>\n";
    *rParamsFile << "\t\t\t<Threshold>" << mThreshold << "</Threshold>\n";

    AbstractForce<DIM, DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CircadianRandomForce<1>;
template class CircadianRandomForce<2>;
template class CircadianRandomForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircadianRandomForce)
