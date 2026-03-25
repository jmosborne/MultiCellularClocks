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

#ifndef TEST2DNODEBASEDNINITIALCELLS_HPP_
#define TEST2DNODEBASEDNINITIALCELLS_HPP_

#include <cmath>
#include <memory>
#include <string>

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CellLabelWriter.hpp"
#include "CircadianRandomForce.hpp"
#include "CircadianRhythmModifier.hpp"
#include "NonCyclingCellCountWriter.hpp"
#include "NonCyclingCellProperty.hpp"
#include "RepulsionForce.hpp"
#include "SphereBasedBoundaryCondition.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "VolumeTrackingModifier.hpp"

// PETSc must be initialized to solve linear algebra problems in Chaste.
#include "FakePetscSetup.hpp"

/**
 * 2D node-based simulation starting from n separate initial cells.
 */
class Test2dNodeBasedNInitialCells : public AbstractCellBasedTestSuite
{
public:
    void TestNodeBased2dFromNSeparateInitialCells()
    {
        EXIT_IF_PARALLEL;

        // Number of initial cells (change this to run with a different n).
        const unsigned n_initial_cells = 64u; //64u;

        // Place cells on a square lattice with spacing >= 1.0 so they start separated.
        const unsigned side = static_cast<unsigned>(std::ceil(std::sqrt(static_cast<double>(n_initial_cells))));
        const double spacing = 1.0;

        for (unsigned sim_index = 0u; sim_index < 10u; ++sim_index)
        {
            std::cout << " Run number " << sim_index << "... \n" << std::flush;

            // Reseed the random number generator
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(sim_index);

            std::vector<std::unique_ptr<Node<2> > > node_storage;
            node_storage.reserve(n_initial_cells);

            std::vector<Node<2>*> nodes;
            nodes.reserve(n_initial_cells);

            for (unsigned i = 0u; i < n_initial_cells; ++i)
            {
                const unsigned row = i / side;
                const unsigned col = i % side;
                std::vector<double> location(2);
                location[0] = spacing * (static_cast<double>(col) - 0.5 * static_cast<double>(side - 1u));
                location[1] = spacing * (static_cast<double>(row) - 0.5 * static_cast<double>(side - 1u));

                node_storage.emplace_back(std::make_unique<Node<2> >(i, location, false));
                nodes.push_back(node_storage.back().get());
            }

            NodesOnlyMesh<2> mesh;
            const double max_interaction_radius = 1.5;
            mesh.ConstructNodesWithoutMesh(nodes, max_interaction_radius);

            std::vector<CellPtr> cells;
            CellsGenerator<BernoulliTrialWithContactInhibitionCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, n_initial_cells);

            TS_ASSERT_EQUALS(cells.size(), n_initial_cells);
            TS_ASSERT_EQUALS(n_initial_cells % 2u, 0u);

            for (auto& p_cell : cells)
            {
              auto* p_cycle_model = static_cast<BernoulliTrialWithContactInhibitionCellCycleModel*>(p_cell->GetCellCycleModel());
              p_cycle_model->SetQuiescentVolumeFraction(0.75);
              p_cycle_model->SetEquilibriumVolume(1.0);
              p_cycle_model->SetDivisionProbability(0.02);
              p_cycle_model->SetMinimumDivisionAge(4.0);
            }

            // Randomly mark exactly 50% of cells as non-cycling (without
            // replacement). These cells will not divide and their circadian
            // rhythm will remain frozen at its t=0 value.
            std::vector<unsigned> cell_indices(cells.size());
            for (unsigned i = 0u; i < cell_indices.size(); ++i)
            {
                cell_indices[i] = i;
            }

            // Fisher-Yates shuffle driven by Chaste RNG
            for (unsigned i = 0u; i < cell_indices.size(); ++i)
            {
                unsigned remaining = static_cast<unsigned>(cell_indices.size() - i);
                unsigned j = i + static_cast<unsigned>(std::floor(p_gen->ranf() * remaining));
                if (j >= cell_indices.size())
                {
                    j = static_cast<unsigned>(cell_indices.size() - 1u);
                }
                std::swap(cell_indices[i], cell_indices[j]);
            }

            unsigned num_non_cycling = static_cast<unsigned>(cells.size() / 2u);
            for (unsigned k = 0u; k < num_non_cycling; ++k)
            {
                cells[cell_indices[k]]->AddCellProperty(boost::make_shared<NonCyclingCellProperty>());
            }

            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.AddCellWriter<CellLabelWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
            cell_population.AddPopulationWriter<NonCyclingCellCountWriter>();
            
            cell_population.SetDampingConstantNormal(1.0);
            cell_population.SetDampingConstantMutant(0.1);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("MultiCellularClocks/NodeBased2dNInitialCells/" + std::to_string(sim_index));
            simulator.SetEndTime(10*24.0);
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(50);

            auto p_force = boost::make_shared<RepulsionForce<2> >();
            p_force->SetMeinekeSpringStiffness(1.0);
            simulator.AddForce(p_force);

            auto p_circadian_modifier = boost::make_shared<CircadianRhythmModifier<2> >();
            p_circadian_modifier->SetCircadianPeriod(24.0);
            p_circadian_modifier->SetPhaseShift(0.0);
            p_circadian_modifier->SetMutationThreshold(0.0);
            simulator.AddSimulationModifier(p_circadian_modifier);

            auto p_volume_modifier = boost::make_shared<VolumeTrackingModifier<2> >();
            simulator.AddSimulationModifier(p_volume_modifier);


            // auto p_circadian_random_force = boost::make_shared<CircadianRandomForce<2> >();
            // p_circadian_random_force->SetAmplitude(0.8);
            // simulator.AddForce(p_circadian_random_force);

            // Circular boundary of radius 10 centred at the origin.
            // c_vector<double, 2> centre = zero_vector<double>(2);
            // const double radius = 5.0;
            // auto p_circle_bc = boost::make_shared<SphereBasedBoundaryCondition<2> >(&cell_population, centre, radius);
            // simulator.AddCellPopulationBoundaryCondition(p_circle_bc);
            
            simulator.Solve();

            // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

        }
    }
};

#endif // TEST2DNODEBASEDNINITIALCELLS_HPP_
