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

#include "slab_tip_trench_position.h"
#include <aspect/geometry_model/box.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    SlabTipTrenchPosition<dim>::execute (TableHandler &statistics)
    {
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(0).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> compositional_values_slab(n_q_points);
      std::vector<double> compositional_values_trench(n_q_points);
      std::vector<Point<dim> > position_values(n_q_points);

      // The most right possible trench position
      double local_trench_position = geometry->get_extents()[0]; 
 
      // The minimal depth of the slab
      double local_slab_depth      = 0.0;

      // The maximal depth of the slab
      const double max_depth = geometry->maximal_depth();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {

            fe_values.reinit (cell);

            // Retrieve position of quadrature points
            position_values = fe_values.get_quadrature_points();

            // Retrieve compositional values for the slab and trench
            fe_values[this->introspection().extractors.compositional_fields[slab]].get_function_values
            (this->get_solution(), compositional_values_slab);
            fe_values[this->introspection().extractors.compositional_fields[trench]].get_function_values
            (this->get_solution(), compositional_values_trench);

            // The vertical distance between quadrature points within an element
            const double delta = position_values[3][1] - position_values[0][1];

            // Obtain current slab tip depth and most left point of trench in a loop over the quadrature points
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                // Depth of quadrature point
                const double depth = geometry->depth(position_values[q]);

                // The compositional field go from 0 to 1, so 0.5 is taken to represent their boundary
                if (compositional_values_slab[q] >= 0.5 && depth > local_slab_depth)
                    local_slab_depth = depth;

                // Retrieve the trench position along the surface, which we take as the original depth 
                // of the sticky air of 50 km
                if (compositional_values_trench[q] >= 0.5 && depth <= (50000.0+delta) && depth >= (50000.0-delta) && position_values[q][0] < local_trench_position)
                    local_trench_position = position_values[q][0];
              }
          }


      // Compute the maximum slab tip depth over all processors
      const  double global_slab_depth =
        Utilities::MPI::max (local_slab_depth, this->get_mpi_communicator());

      AssertThrow(global_slab_depth<=max_depth+1000.,ExcMessage("Computed slab tip depth bigger than domain depth."));
   
      // Compute the minimum trench position over all processors
      // by taking the maximum value of the negative trench position
      local_trench_position = -1.0 * local_trench_position;
      const  double global_trench_position =
        Utilities::MPI::max (local_trench_position, this->get_mpi_communicator());

      // Add output to the statistics file
      statistics.add_value ("Slab tip depth [km]", global_slab_depth/1000.0);
      statistics.add_value ("Trench position [km]", -global_trench_position/1000.0);

      // Also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Slab tip depth [km]" ,
                                  "Trench position [km]"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      // Add output to the log file
      std::ostringstream output;
      output.precision(3);
      output << global_slab_depth/1000.0
             << " km "
             << -global_trench_position/1000.0
             << " km";
      return std::pair<std::string, std::string> ("Slab tip depth and trench position [km]:",
                                                  output.str());
    }

    template <int dim>
    void
    SlabTipTrenchPosition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Subduction statistics");
        {
          prm.declare_entry("Compositional field number of slab", "1",
                            Patterns::Integer (0),
                            "The number of the compositional field that represents the slab tip. Units: -.");
          prm.declare_entry("Compositional field number of trench", "2",
                            Patterns::Integer (0),
                            "The number of the compositional field that represents the trench. Units: -.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SlabTipTrenchPosition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Subduction statistics");
        {
          slab   = prm.get_integer("Compositional field number of slab");
          trench = prm.get_integer("Compositional field number of trench");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SlabTipTrenchPosition,
                                  "subduction statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "subduction plate. "
                                  "It computest the deepest point of the compositional field that is identified as "
                                  "representing the tip of the slab in the input file as well as the most left "
                                  "horizontal coordinate of the field identified as representing the start of the trench.")
  }
}
