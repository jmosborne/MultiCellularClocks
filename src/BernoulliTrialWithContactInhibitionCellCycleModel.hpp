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

#ifndef BERNOULLITRIALWITHCONTACTINHIBITIONCELLCYCLEMODEL_HPP_
#define BERNOULLITRIALWITHCONTACTINHIBITIONCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * Bernoulli-trial cell-cycle model with contact inhibition.
 *
 * All non-differentiated cells divide stochastically with a single per-hour
 * probability, once they are older than a single minimum division age.
 *
 * If a cell's current volume is below
 * mQuiescentVolumeFraction * mEquilibriumVolume,
 * it is contact inhibited and cannot divide.
 */
class BernoulliTrialWithContactInhibitionCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mQuiescentVolumeFraction;
        archive & mEquilibriumVolume;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;
        archive & mDivisionProbability;
        archive & mMinimumDivisionAge;
    }

protected:
    double mQuiescentVolumeFraction;
    double mEquilibriumVolume;

    double mCurrentQuiescentOnsetTime;
    double mCurrentQuiescentDuration;

    double mDivisionProbability;
    double mMinimumDivisionAge;

    BernoulliTrialWithContactInhibitionCellCycleModel(const BernoulliTrialWithContactInhibitionCellCycleModel& rModel)
        : AbstractCellCycleModel(rModel),
          mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
          mEquilibriumVolume(rModel.mEquilibriumVolume),
          mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
          mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration),
            mDivisionProbability(rModel.mDivisionProbability),
            mMinimumDivisionAge(rModel.mMinimumDivisionAge)
    {
    }

public:
    BernoulliTrialWithContactInhibitionCellCycleModel()
        : AbstractCellCycleModel(),
          mQuiescentVolumeFraction(DOUBLE_UNSET),
          mEquilibriumVolume(DOUBLE_UNSET),
          mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
          mCurrentQuiescentDuration(0.0),
            mDivisionProbability(0.1),
            mMinimumDivisionAge(1.0)
    {
    }

    bool ReadyToDivide() override
    {
        assert(mpCell != nullptr);

        if (!mReadyToDivide)
        {
            // Differentiated cells do not divide.
            if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
            {
                return false;
            }

            if ((mQuiescentVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
            {
                EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
            }

            // Contact inhibition check.
            double cell_volume = mpCell->GetCellData()->GetItem("volume");
            double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

            if (cell_volume < quiescent_volume)
            {
                mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;

                // Label contact-inhibited cells.
                boost::shared_ptr<AbstractCellProperty> p_label =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
                mpCell->AddCellProperty(p_label);

                return false;
            }

            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();

            // Cell is not contact inhibited, so remove label if present.
            if (mpCell->template HasCellProperty<CellLabel>())
            {
                mpCell->template RemoveCellProperty<CellLabel>();
            }

            double dt = SimulationTime::Instance()->GetTimeStep();
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

            if (GetAge() > mMinimumDivisionAge)
            {
                if (p_gen->ranf() < mDivisionProbability*dt)
                {
                    mReadyToDivide = true;
                }
            }
        }

        return mReadyToDivide;
    }

    AbstractCellCycleModel* CreateCellCycleModel() override
    {
        return new BernoulliTrialWithContactInhibitionCellCycleModel(*this);
    }

    void SetQuiescentVolumeFraction(double quiescentVolumeFraction)
    {
        mQuiescentVolumeFraction = quiescentVolumeFraction;
    }

    double GetQuiescentVolumeFraction() const
    {
        return mQuiescentVolumeFraction;
    }

    void SetEquilibriumVolume(double equilibriumVolume)
    {
        mEquilibriumVolume = equilibriumVolume;
    }

    double GetEquilibriumVolume() const
    {
        return mEquilibriumVolume;
    }

    double GetCurrentQuiescentDuration() const
    {
        return mCurrentQuiescentDuration;
    }

    double GetCurrentQuiescentOnsetTime() const
    {
        return mCurrentQuiescentOnsetTime;
    }

    void SetDivisionProbability(double divisionProbability)
    {
        mDivisionProbability = divisionProbability;
    }

    double GetDivisionProbability()
    {
        return mDivisionProbability;
    }

    void SetMinimumDivisionAge(double minimumDivisionAge)
    {
        mMinimumDivisionAge = minimumDivisionAge;
    }

    double GetMinimumDivisionAge()
    {
        return mMinimumDivisionAge;
    }

    double GetAverageTransitCellCycleTime()
    {
        return mMinimumDivisionAge + 1.0/mDivisionProbability;
    }

    double GetAverageStemCellCycleTime()
    {
        return mMinimumDivisionAge + 1.0/mDivisionProbability;
    }

    void OutputCellCycleModelParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
        *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";
        *rParamsFile << "\t\t\t<DivisionProbability>" << mDivisionProbability << "</DivisionProbability>\n";
        *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

        AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    }
};

#endif // BERNOULLITRIALWITHCONTACTINHIBITIONCELLCYCLEMODEL_HPP_
