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

#ifndef SPHEREBASEDBOUNDARYCONDITION_HPP_
#define SPHEREBASEDBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "ChasteSerialization.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A spherical boundary condition that keeps cells inside (or on) a sphere.
 *
 * In 2D, this corresponds to keeping cells inside a circle.
 */
template<unsigned DIM>
class SphereBasedBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:
    /** Centre of the sphere. */
    c_vector<double, DIM> mCentre;

    /** Radius of the sphere. */
    double mRadius;

    /** Numerical tolerance used by VerifyBoundaryCondition. */
    double mTolerance;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mRadius;
        archive & mTolerance;
    }

public:
    SphereBasedBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                 c_vector<double, DIM> centre,
                                 double radius,
                                 double tolerance=1e-8)
        : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
          mCentre(centre),
          mRadius(radius),
          mTolerance(tolerance)
    {
        assert(mRadius > 0.0);
        assert(mTolerance > 0.0);

        if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == nullptr)
        {
            EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
        }
    }

    virtual ~SphereBasedBoundaryCondition()
    {
    }

    const c_vector<double, DIM>& rGetCentre() const
    {
        return mCentre;
    }

    double GetRadius() const
    {
        return mRadius;
    }

    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations) override
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            c_vector<double, DIM> offset = cell_location - mCentre;
            double distance = norm_2(offset);

            if (distance > mRadius)
            {
                c_vector<double, DIM> new_location = mCentre + (mRadius / distance) * offset;
                unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                p_node->rGetModifiableLocation() = new_location;
            }
        }
    }

    bool VerifyBoundaryCondition() override
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double distance = norm_2(cell_location - mCentre);
            if (distance > mRadius + mTolerance)
            {
                return false;
            }
        }
        return true;
    }

    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<Centre>";
        for (unsigned i=0; i != DIM-1U; i++)
        {
            *rParamsFile << mCentre[i] << ",";
        }
        *rParamsFile << mCentre[DIM-1] << "</Centre>\n";
        *rParamsFile << "\t\t\t<Radius>" << mRadius << "</Radius>\n";
        *rParamsFile << "\t\t\t<Tolerance>" << mTolerance << "</Tolerance>\n";

        AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SphereBasedBoundaryCondition)

#endif // SPHEREBASEDBOUNDARYCONDITION_HPP_
