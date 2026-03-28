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

#include "CiradianRepulsionForce.hpp"

#include <cassert>
#include <cmath>

template<unsigned DIM>
CiradianRepulsionForce<DIM>::CiradianRepulsionForce()
    : AbstractForce<DIM, DIM>(),
    mRepulsionParameter(15.0)
{
}

template<unsigned DIM>
CiradianRepulsionForce<DIM>::~CiradianRepulsionForce()
{
}

template<unsigned DIM>
double CiradianRepulsionForce<DIM>::GetRepulsionParameter() const
{
    return mRepulsionParameter;
}

template<unsigned DIM>
void CiradianRepulsionForce<DIM>::SetRepulsionParameter(double repulsionParameter)
{
    assert(repulsionParameter > 0.0);
    mRepulsionParameter = repulsionParameter;
}

template<unsigned DIM>
void CiradianRepulsionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    auto* p_node_based_population = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (p_node_based_population == nullptr)
    {
        EXCEPTION("CiradianRepulsionForce is to be used with a NodeBasedCellPopulation only");
    }

    const std::vector<std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = p_node_based_population->rGetNodePairs();

    for (const auto& [p_node_a, p_node_b] : r_node_pairs)
    {
        const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
        const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

        const double node_a_radius = p_node_a->GetRadius();
        const double node_b_radius = p_node_b->GetRadius();
        const double rest_length = node_a_radius + node_b_radius;

        c_vector<double, DIM> vector_a_to_b = p_node_based_population->rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);
        const double distance_between_nodes = norm_2(vector_a_to_b);
        assert(distance_between_nodes > 0.0);

        if (distance_between_nodes < rest_length)
        {
            c_vector<double, DIM> unit_difference = vector_a_to_b / distance_between_nodes;
            const double overlap = distance_between_nodes - rest_length;

            // Same repulsive branch as GeneralisedLinearSpringForce for overlap <= 0.
            assert(overlap > -rest_length);
            c_vector<double, DIM> force = mRepulsionParameter * unit_difference * rest_length * std::log(1.0 + overlap / rest_length);

            for (unsigned j = 0; j < DIM; ++j)
            {
                assert(!std::isnan(force[j]));
            }

            // Read per-cell drag coefficients from CellData and weight force
            // contributions by inverse drag (higher drag => smaller force contribution).
            double drag_a = 1.0;
            double drag_b = 1.0;

            CellPtr p_cell_a = p_node_based_population->GetCellUsingLocationIndex(p_node_a->GetIndex());
            CellPtr p_cell_b = p_node_based_population->GetCellUsingLocationIndex(p_node_b->GetIndex());

            if (p_cell_a->GetCellData()->HasItem("damping_coefficient"))
            {
                drag_a = p_cell_a->GetCellData()->GetItem("damping_coefficient");
            }

            if (p_cell_b->GetCellData()->HasItem("damping_coefficient"))
            {
                drag_b = p_cell_b->GetCellData()->GetItem("damping_coefficient");
            }

            p_node_a->AddAppliedForceContribution((1.0 / drag_a) * force);
            p_node_b->AddAppliedForceContribution((-1.0 / drag_b) * force);
        }
    }
}

template<unsigned DIM>
void CiradianRepulsionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RepulsionParameter>" << mRepulsionParameter << "</RepulsionParameter>\n";

    AbstractForce<DIM, DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CiradianRepulsionForce<1>;
template class CiradianRepulsionForce<2>;
template class CiradianRepulsionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CiradianRepulsionForce)
