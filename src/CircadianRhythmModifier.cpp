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

#include "CircadianRhythmModifier.hpp"

#include <cassert>
#include <cmath>

#include "NonCyclingCellProperty.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
CircadianRhythmModifier<DIM>::CircadianRhythmModifier()
    : AbstractCellBasedSimulationModifier<DIM, DIM>(),
      mCircadianPeriod(24.0),
      mPhaseShift(0.0)
{
}

template<unsigned DIM>
CircadianRhythmModifier<DIM>::~CircadianRhythmModifier()
{
}

template<unsigned DIM>
double CircadianRhythmModifier<DIM>::GetCircadianPeriod() const
{
    return mCircadianPeriod;
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::SetCircadianPeriod(double circadianPeriod)
{
    assert(circadianPeriod > 0.0);
    mCircadianPeriod = circadianPeriod;
}

template<unsigned DIM>
double CircadianRhythmModifier<DIM>::GetPhaseShift() const
{
    return mPhaseShift;
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::SetPhaseShift(double phaseShift)
{
    mPhaseShift = phaseShift;
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                                              std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    const double time = SimulationTime::Instance()->GetTime();
    const double circadian_cycle = std::sin(2.0 * 3.14159265358979323846 * time / 24.0);

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator pCell = rCellPopulation.Begin();
         pCell != rCellPopulation.End();
         ++pCell)
    {
        // Non-cycling cells have a frozen circadian rhythm: their CellData
        // items are written once on the first call (SetupSolve at t=0) and
        // are not updated on subsequent time steps.
        if (pCell->template HasCellProperty<NonCyclingCellProperty>()
            && pCell->GetCellData()->HasItem("circadian_cycle"))
        {
            continue;
        }

        pCell->GetCellData()->SetItem("circadian_cycle", circadian_cycle);
    }
}

template<unsigned DIM>
void CircadianRhythmModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CircadianPeriod>" << mCircadianPeriod << "</CircadianPeriod>\n";
    *rParamsFile << "\t\t\t<PhaseShift>" << mPhaseShift << "</PhaseShift>\n";

    AbstractCellBasedSimulationModifier<DIM, DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CircadianRhythmModifier<1>;
template class CircadianRhythmModifier<2>;
template class CircadianRhythmModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircadianRhythmModifier)
