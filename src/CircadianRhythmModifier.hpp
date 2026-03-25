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

#ifndef CIRCADIANRHYTHMMODIFIER_HPP_
#define CIRCADIANRHYTHMMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * Assigns a circadian rhythm variable to each cell.
 *
 * The following CellData item is updated every time step:
 *  - "circadian_cycle" = sin(2*pi*time/24)
 */
template<unsigned DIM>
class CircadianRhythmModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
        archive & mCircadianPeriod;
        archive & mPhaseShift;
    }

    /** Circadian period in simulation time units (default: 24). */
    double mCircadianPeriod;

    /** Optional phase shift in simulation time units (default: 0). */
    double mPhaseShift;

public:

    /** Constructor. */
    CircadianRhythmModifier();

    /** Destructor. */
    virtual ~CircadianRhythmModifier();

    /**
     * @return circadian period.
     */
    double GetCircadianPeriod() const;

    /**
     * Set circadian period.
     *
     * @param circadianPeriod period (must be > 0)
     */
    void SetCircadianPeriod(double circadianPeriod);

    /**
     * @return phase shift.
     */
    double GetPhaseShift() const;

    /**
     * Set phase shift.
     *
     * @param phaseShift phase shift in time units
     */
    void SetPhaseShift(double phaseShift);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override;

    /**
     * Overridden SetupSolve() method.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory output directory relative to Chaste output path
     */
    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory) override;

    /**
     * Updates circadian data for all cells.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircadianRhythmModifier)

#endif // CIRCADIANRHYTHMMODIFIER_HPP_
