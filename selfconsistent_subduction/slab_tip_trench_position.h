/*
  Copyright (C) 2011-2017 by the authors of the ASPECT code.

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

#ifndef __aspect__postprocess_slab_tip_trench_position_h
#define __aspect__postprocess_slab_tip_trench_position_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the deepest point of the composition
     * that represents the subducting plate (slab) and it's most left horizontal position at the
     * surface, i.e. the position of the trench.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class SlabTipTrenchPosition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
       /**
        * Evaluate the solution for the deepest and most left position of the slab.
        */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

      /**
       * Declare the parameters this class takes through input files.
       */
       static
       void
       declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
       virtual
       void
       parse_parameters (ParameterHandler &prm);

      private:
      /**
       * The compositional field id of the field that represents the subducting plate.
       */
       int slab;

      /**
       * The compositional field id of the field that represents the subducting plate at the surface.
       */
       int trench;

    };
  }
}


#endif
