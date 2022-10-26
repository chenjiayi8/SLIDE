/*
 * Module_base_p.hpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "Module.hpp"
#include "../types/data_storage/cell_data.hpp"
#include "../utility/utility.hpp"

#include <string>

namespace slide {
class Module_p : public Module
{
protected:
#if TIMING //!< std::clock_t tstart;
  TimingData_Module_p timeData{};
#endif
public:
  Module_p() : Module("moduleP") {} //!< note this constructor should never be used. It can't determine which coolsystem to use
  Module_p(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype)
    : Module(ID_, Ti, print, pari, Ncells, coolControl, cooltype) {}

  //!< functions from Module_base
  //!< the voltage limits are the most constraining limits of all cells
  //!< ie the highest Vmin of SUs is the Vmin of the module
  double Vmin() override { return transform_max(SUs, free::get_Vmin<SU_t>); }
  double VMIN() override { return transform_max(SUs, free::get_VMIN<SU_t>); }

  double Vmax() override { return transform_min(SUs, free::get_Vmax<SU_t>); }
  double VMAX() override { return transform_min(SUs, free::get_VMAX<SU_t>); }

  double I() override { return transform_sum(SUs, free::get_I<SU_t>); }     //!< the current is the sum  of the current of each cell. Returns 0 if empty.
  double Cap() override { return transform_sum(SUs, free::get_Cap<SU_t>); } //!< module capacity (sum of cells)

  // #TODO -> getOCV print parameter is inactive here...
  double getOCV(bool print = true) override { return transform_mean(SUs, free::get_OCV<SU_t>); }
  double getRtot() override;
  double V(bool print = true) override;      //!< module voltage (sum of cells), print is an optional argument
  double getVi(size_t i, bool print = true); //!< get the voltage of SU[i] while accounting for the contact resistance

  Status setI_iterative(double Inew, bool checkV = true, bool print = true);      //!< set a module current iteratively
  Status setCurrent(double Inew, bool checkV = true, bool print = true) override; //!< set a module current
  Status redistributeCurrent(bool checkV = true, bool print = true);              //!< redistribute the total module current to the different cells

  bool validSUs(SUs_span_t c, bool print = true) override; //!< check if the cells in this array are valid for this module

  void timeStep_CC(double dt, int steps = 1) override;

  Module_p *copy() override;
  TimingData_Module_p getTimings();
  void setTimings(TimingData_Module_p td);
};
} // namespace slide
