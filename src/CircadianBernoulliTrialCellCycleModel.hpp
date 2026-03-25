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

#ifndef CIRCADIANBERNOULLITRIALCELLCYCLEMODEL_HPP_
#define CIRCADIANBERNOULLITRIALCELLCYCLEMODEL_HPP_

#include <cmath>

#include "BernoulliTrialCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

/**
 * Bernoulli trial cell-cycle model in which the division probability is
 * modulated by the circadian cycle angle stored in CellData.
 *
 * The instantaneous per-hour division probability is
 *
 * p_div = p0 * (1 + A * cos(theta + phi))
 *
 * where
 * - p0 is the baseline probability (mDivisionProbability),
 * - A is mCircadianModulationAmplitude,
 * - theta is the cell-data item "circadian_cycle_angle",
 * - phi is mCircadianAngleOffset.
 */
class CircadianBernoulliTrialCellCycleModel : public BernoulliTrialCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing.
     *
     * @param archive The boost archive.
     * @param version The current version of this class.
     */
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<BernoulliTrialCellCycleModel>(*this);
        archive & mCircadianModulationAmplitude;
        archive & mCircadianAngleOffset;
    }

    /** Modulation amplitude A. */
    double mCircadianModulationAmplitude;

    /** Optional angle offset phi in radians. */
    double mCircadianAngleOffset;

protected:

    /**
     * Copy constructor used by CreateCellCycleModel().
     *
     * @param rModel model to copy
     */
    CircadianBernoulliTrialCellCycleModel(const CircadianBernoulliTrialCellCycleModel& rModel)
        : BernoulliTrialCellCycleModel(rModel),
          mCircadianModulationAmplitude(rModel.mCircadianModulationAmplitude),
          mCircadianAngleOffset(rModel.mCircadianAngleOffset)
    {
    }

public:

    /** Constructor. */
    CircadianBernoulliTrialCellCycleModel()
        : BernoulliTrialCellCycleModel(),
          mCircadianModulationAmplitude(0.5),
          mCircadianAngleOffset(0.0)
    {
    }

    /**
     * Overridden ReadyToDivide() method with circadian modulation.
     *
     * @return whether the cell is ready to divide
     */
    bool ReadyToDivide() override
    {
        assert(mpCell != nullptr);

        if (!mReadyToDivide)
        {
            if (GetAge() > mMinimumDivisionAge)
            {
                if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
                {
                    double angle = 0.0;
                    if (mpCell->GetCellData()->HasItem("circadian_cycle_angle"))
                    {
                        angle = mpCell->GetCellData()->GetItem("circadian_cycle_angle");
                    }

                    const double modulation = 1.0 + mCircadianModulationAmplitude*std::cos(angle + mCircadianAngleOffset);
                    double effective_division_probability = mDivisionProbability*modulation;
                    if (effective_division_probability < 0.0)
                    {
                        effective_division_probability = 0.0;
                    }

                    const double dt = SimulationTime::Instance()->GetTimeStep();
                    if (RandomNumberGenerator::Instance()->ranf() < effective_division_probability*dt)
                    {
                        mReadyToDivide = true;
                    }
                }
            }
        }

        return mReadyToDivide;
    }

    /**
     * Overridden builder method to create new instances of the model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel() override
    {
        return new CircadianBernoulliTrialCellCycleModel(*this);
    }

    /** @return mCircadianModulationAmplitude */
    double GetCircadianModulationAmplitude() const
    {
        return mCircadianModulationAmplitude;
    }

    /**
     * Set mCircadianModulationAmplitude.
     *
     * @param amplitude modulation amplitude
     */
    void SetCircadianModulationAmplitude(double amplitude)
    {
        mCircadianModulationAmplitude = amplitude;
    }

    /** @return mCircadianAngleOffset */
    double GetCircadianAngleOffset() const
    {
        return mCircadianAngleOffset;
    }

    /**
     * Set mCircadianAngleOffset.
     *
     * @param angleOffset angle offset in radians
     */
    void SetCircadianAngleOffset(double angleOffset)
    {
        mCircadianAngleOffset = angleOffset;
    }

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile stream to output parameter values
     */
    void OutputCellCycleModelParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<CircadianModulationAmplitude>" << mCircadianModulationAmplitude
                     << "</CircadianModulationAmplitude>\n";
        *rParamsFile << "\t\t\t<CircadianAngleOffset>" << mCircadianAngleOffset
                     << "</CircadianAngleOffset>\n";

        BernoulliTrialCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CircadianBernoulliTrialCellCycleModel)

#endif // CIRCADIANBERNOULLITRIALCELLCYCLEMODEL_HPP_
