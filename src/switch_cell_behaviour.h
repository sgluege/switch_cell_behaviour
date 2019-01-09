// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef SWITCH_CELL_BEHAVIOUR_H_
#define SWITCH_CELL_BEHAVIOUR_H_

#include "biodynamo.h"

// set numer of simulation steps
const int simulation_steps = 500; // Time between two simulation steps equals: 0.01hours (default)

// set amount of intracellular_substance_quantity for precursor cells
const double intracellular_substance_quantity = 100;
const double counter_sub_quant_threshold = 75; // threshold of intracellular substance ... when to switch biomodules
const double default_cell_diameter = 6; // default diameter of new cells

namespace bdm {

// 1. Define general cell Module
    struct GeneralModule : public BaseBiologyModule {
        GeneralModule() : BaseBiologyModule(gAllEventIds) {}

        /// Empty default event constructor, because GrowthModule does not have state.
        template <typename TEvent, typename TBm>
        GeneralModule(const TEvent& event, TBm* other, uint64_t new_oid = 0)
                : BaseBiologyModule(event, other, new_oid) {}

        template <typename T, typename TSimulation = Simulation<>>
        void Run(T* cell) {

            // read current cell type
            int current_cell_type = cell->GetCellType();

            // run cell behaviour depending on cell type
            if (current_cell_type==1){
                type1behaviour(cell);
            } else if (current_cell_type==2){
                type2behaviour(); // yet ... do nothing
            } else {
                std::cout << "no behaviour defined for cell of type: " << current_cell_type << std::endl;
            }
        }

        // define behaviour of cell type S1
        template <typename T, typename TSimulation = Simulation<>>
        void type1behaviour(T* cell){

            // change counter_sub_quant_ (intracellular substance quantity) ... reduce by 1% per timestep
            double currentCounterSubQuant = cell->GetCounterSubQuant();
            currentCounterSubQuant = currentCounterSubQuant * (0.99);
            cell->SetCounterSubQuant(currentCounterSubQuant);

            // cell differentiation
            runCellCycleDiffStep(cell); // grow cell

            // intracellular substance quantity is below threshold
            // so switch to cell type 2
            if (currentCounterSubQuant < counter_sub_quant_threshold) {
                cell->SetCellType(2);
                cell->SetCellColor(1);
                cell->SetCounterSubQuant(0);
                std::cout << "switch to cell type 2" << std::endl;
            }

        }

        // define behaviour of cell type S2
        void type2behaviour(){
            // yet ... do nothing
        }

        // run cell differentiation step (used in type1)
        template <typename T, typename TSimulation = Simulation<>>
        void runCellCycleDiffStep(T* cell){
            if (cell->GetDiameter() < cell->GetMaxDiam()) { // max_diam not reached
                cell->SetDiameter(cell->GetDiameter() + 0.2); // grow cell diameter by adding some number
            } else { // max diam reached -> divide cell
                auto &&daughter = cell->Divide();
                daughter->SetDiameter(default_cell_diameter);  // init daughter with default diameter
            }
        }

    private:
        BDM_CLASS_DEF_NV(GeneralModule, 1);
    };


// Define new cell (BDM_SIM_OBJECT) which extends Cell by adding extra data
BDM_SIM_OBJECT(MyCell, Cell) {  // our object extends the Cell object
// create the header with our new data member
BDM_SIM_OBJECT_HEADER(MyCell, Cell, 1, cell_type_, cell_color_, counter_sub_quant_, max_diam_);

public:
    MyCellExt() {}
    explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {
        // init cell with type 1, color 0 and max_diam = 10
        cell_type_[kIdx] = 1;
        cell_color_[kIdx] = 0;
        counter_sub_quant_[kIdx] = intracellular_substance_quantity; // init with type 1 and specific intracellular_substance_quantity
        max_diam_[kIdx] = 10;
    }

    /// If MyCell divides, daughter 2 copies the data members from the mother
    template <typename TMother>
    MyCellExt(const CellDivisionEvent& event, TMother* mother)
            : Base(event, mother) {
        cell_type_[kIdx] = mother->cell_type_[mother->kIdx];
        cell_color_[kIdx] = mother->cell_color_[mother->kIdx];
        counter_sub_quant_[kIdx] = mother->counter_sub_quant_[mother->kIdx]; // intracellular substance quantity
        max_diam_[kIdx] = mother->max_diam_[mother->kIdx];
    }

    /// If a cell divides, daughter keeps the same state from its mother.
    template <typename TDaughter>
    void EventHandler(const CellDivisionEvent& event, TDaughter* daughter) {
        Base::EventHandler(event, daughter);
    }

    // getter and setter for our new data member
    void SetCellType(int cell_type) { cell_type_[kIdx] = cell_type; }
    int GetCellType() const { return cell_type_[kIdx]; }

    void SetCellColor(int cell_color) { cell_color_[kIdx] = cell_color; }
    int GetCellColor() const { return cell_color_[kIdx]; }

    void SetCounterSubQuant(double counter_sub_quant) { counter_sub_quant_[kIdx] = counter_sub_quant; }
    double GetCounterSubQuant() const { return counter_sub_quant_[kIdx]; }

    void SetMaxDiam(double max_diam) { max_diam_[kIdx] = max_diam; }
    double GetMaxDiam() const { return max_diam_[kIdx]; }

private:
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    vec<int> cell_type_;
    vec<int> cell_color_;
    vec<double> counter_sub_quant_;  // intracellular substance quantity
    vec<double> max_diam_;
};


// Define compile time parameter
BDM_CTPARAM() {
    BDM_CTPARAM_HEADER();
    using SimObjectTypes = CTList<MyCell>;  // use MyCell object
    // Override default BiologyModules for Cell
    BDM_CTPARAM_FOR(bdm, MyCell) { using BiologyModules = CTList<GeneralModule>; };

};

inline int Simulate(int argc, const char** argv) {
    // set parameters of the simulation
    auto set_param = [](auto* param) {
        param->bound_space_ = true;
        param->min_bound_ = -50;
        param->max_bound_ = 50;  // plate of 100*100
        param->run_mechanical_interactions_ = true;
    };

    Simulation<> simulation(argc, argv, set_param);

    // Define initial model - in this example: single cell at origin
    auto* rm = simulation.GetResourceManager();

    // create a structure to contain cells
    auto* cells = rm->template Get<MyCell>();
    // allocate the correct number of cell in our cells structure before
    cells->reserve(1);

    MyCell cell({0, 0, 10});
    cell.SetDiameter(default_cell_diameter);
    cell.AddBiologyModule(GeneralModule());
    cells->push_back(cell);  // put the created cell in cells structure
    cells->Commit();  // commit cells

    // Run simulation
    simulation.GetScheduler()->Simulate(simulation_steps);
    std::cout << "Simulation completed successfully!" << std::endl;
    return 0;
}


}  // namespace bdm

#endif  // SWITCH_CELL_BEHAVIOUR_H_
