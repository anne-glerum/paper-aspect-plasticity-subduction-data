/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__initial_conditions_2Dsubd_temp_h
#define __aspect__initial_conditions_2Dsubd_temp_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes the temperature field according to the plate cooling model
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class SubdTemp : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         *The parameters needed for the plate cooling model
         */ 

        /**
         * The maximum thickness and temperature of an oceanic plate
         * when time goes to infinity
         */   
        double d_max;
        double Tm;

        /**
         * A function object representing the compositional fields
         * that will be used as a reference profile for calculating
         * the thermal diffusivity.
         * The function depends only on depth.
         */
//        std::auto_ptr<Functions::ParsedFunction<dim> > function;
        /**
         * The age of the plates at the trench
         */
        double age_OP_max;
        double age_SP_max;

        /**
         * The  thickness of the oceanic crust
         */
        double d_crust;

        /**
         * The horizontal position of the trench (beginning of WZ).
         */
        double x_trench;

        /**
         * The angle of the slab tip.
         */
        double slab_dip;
    };
  }
}

#endif
