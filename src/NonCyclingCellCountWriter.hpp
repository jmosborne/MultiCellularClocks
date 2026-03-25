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

#ifndef NONCYCLINGCELLCOUNTWRITER_HPP_
#define NONCYCLINGCELLCOUNTWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ChasteSerialization.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NonCyclingCellProperty.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Population writer recording the number of cells with NonCyclingCellProperty
 * and the number of cells without it.
 *
 * Each output line has the form:
 *   [time]\t[num_non_cycling]\t[num_other]
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class NonCyclingCellCountWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

    void VisitAnyPopulation(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned num_non_cycling = 0u;
        unsigned num_other = 0u;

        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
             cell_iter != pCellPopulation->End();
             ++cell_iter)
        {
            if (cell_iter->template HasCellProperty<NonCyclingCellProperty>())
            {
                ++num_non_cycling;
            }
            else
            {
                ++num_other;
            }
        }

        *this->mpOutStream << num_non_cycling << "\t" << num_other;
    }

public:
    NonCyclingCellCountWriter()
        : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("results.non_cycling_cell_counts")
    {
    }

    void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }

    void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }

    void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }

    void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }

    void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }

    void Visit(ImmersedBoundaryCellPopulation<SPACE_DIM>* pCellPopulation) override
    {
        VisitAnyPopulation(pCellPopulation);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NonCyclingCellCountWriter)

#endif // NONCYCLINGCELLCOUNTWRITER_HPP_
