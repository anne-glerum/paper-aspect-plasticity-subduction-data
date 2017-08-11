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


#include "3Dsubd_compo.h"
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/boundary_temperature/box.h>

namespace aspect
{
  namespace CompositionalInitialConditions
  {

    template <int dim>
    double
    SubdCompo<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // this initial condition only makes sense if the geometry is a
      // Box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      Assert (geometry != 0,
              ExcMessage ("This initial condition can only be used if the geometry "
                          "is a box."));

      // this initial condition only makes sense if the boundary temperature model is a
      // Box. verify that it is indeed
      const BoundaryTemperature::Box<dim> *T_boundary
        = dynamic_cast<const BoundaryTemperature::Box<dim>*> (&this->get_boundary_temperature());
      Assert (T_boundary != 0,
              ExcMessage ("This initial condition can only be used if the boundary temperature model "
                          "is a box."));
    
      Assert (dim == 3, ExcMessage ("This initial condition should be used in 3D."));

      //Retrieve the spreading velocities of the plates in seconds if necessary
      //From this we will determine the age of the plate with distance from the ridge
      const double v_spread_OP = (this->convert_output_to_years() ? v_OP / year_in_seconds
                                 : v_OP);
      const double v_spread_SP = (this->convert_output_to_years() ? v_SP / year_in_seconds
                                 : v_SP);
      const double age_AP      = (this->convert_output_to_years() ? t_AP * year_in_seconds
                                 : t_AP);


      const double depth = geometry->depth(position);
      double kap = 1.0e-6;

      //Determine plate age based on distance from domain boundary
      const double age_OP = (position[0]<5000.0) ? 0.0 : (position[0]-5000.0) * (1.0/v_spread_OP);
      const double age_SP = (position[0]>geometry->get_extents()[0]-5000.0)? 0.0 : (geometry->get_extents()[0]-position[0]-5000.0) * (1.0/v_spread_SP);

      // The height of the domain and of the maximum plate thickness
      const double domain_height = geometry->get_extents()[dim-1];
      const double z_max = domain_height-d_max;
      const double domain_width = geometry->get_extents()[0];

      // The slab dip angle in radians
      const double slab_angle = 29.0 * numbers::PI / 180.0;

      // Trench position
      const double trench = 1600000.0;

      // Horizontal (x) position of the slab tip
      const double slab_tip = trench - std::cos(slab_angle) * 206000.0;


      // The SP
      if (position[0] >= trench  && \
          position[dim-2] > plate_width) 
         { 
          // crust
          if (position[dim-1] >= 649600.0 && \
              position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_SP))+20800.0)
             return (n_comp == 4 ) ? 1.0 : 0.0;
          // mantle
          // if we're at the vertical domain wall, don't set composition
          else if (position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_SP)))
             return (n_comp == 2 ) ? (age_SP == 0 ? 0.0 : 1.0) : 0.0;
          else 
          return 0.0;
         }
      // The slab tip
      else if (position[0] >= slab_tip && \
               position[0] < trench &&  \
               position[dim-1] <= domain_height-std::tan(slab_angle) * (trench-position[0]) && \
               position[dim-2] > plate_width)
              {
               //Prescribe T in slab tip by transposing the top of the slab to the surface
               Point<dim> new_position; 
               new_position[0] = position[0];
               new_position[dim-1] = position[dim-1]+std::tan(slab_angle) * (trench-position[0]);
               const double new_depth = geometry->depth(new_position);

               if (new_position[dim-1] >= 649600.0)
                  return (n_comp == 4 ) ? 2.0 : 0.0;
               else if (new_position[dim-1] < 649600.0 && \
                        new_position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_SP)) )
                  return (n_comp == 2) ? 1.0 : 0.0;
               else
                  return 0.0;
              }
      // The OP
      else if (position[dim-1] > domain_height-std::tan(slab_angle) * (trench-position[0]) && \
               position[dim-2] > plate_width)
              {
               if (position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_OP)))
                  return (n_comp == 1 ) ? (age_OP == 0.0 ? 0.0 : 1.0) : 0.0;
               else
                  return 0.0;
              }
      // The AP
      else if (position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_AP)) && position[dim-2] <= plate_width-40000.0)
      {
       return (n_comp == 0) ? 1.0 : 0.0;
      }
      // The WZ
      else if (position[dim-2] > plate_width-40000.0 && position[dim-2] <= plate_width)
      {
       if ((position[0] >= trench && position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_SP))) || \
           (position[0] < trench && position[dim-1] >= domain_height-(2.32*std::sqrt(kap*age_OP))) )
          return (n_comp == 3) ? 1.0 : 0.0;
       else
          return 0.0;
      }      
      else
       return 0.0;


    }

    template <int dim>
    void
    SubdCompo<dim>::declare_parameters (ParameterHandler &prm)
    {
    }


    template <int dim>
    void
    SubdCompo<dim>::parse_parameters (ParameterHandler &prm)
    {
      // we need to get at the number of compositional fields here to
      // initialize the function parser. unfortunately, we can't get it
      // via SimulatorAccess from the simulator itself because at the
      // current point the SimulatorAccess hasn't been initialized
      // yet. so get it from the parameter file directly.
      prm.enter_subsection ("Compositional fields");
      const unsigned int n_compositional_fields = prm.get_integer ("Number of fields");
      prm.leave_subsection ();

      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plate cooling");
        {
          v_OP = prm.get_double ("Spreading velocity overriding plate");
          v_SP = prm.get_double ("Spreading velocity subducting plate");
          plate_width = prm.get_double ("Width of subducting and overriding plate");
          d_max = prm.get_double ("Maximum oceanic plate thickness");
          t_AP = prm.get_double ("Age adjacent plate");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace CompositionalInitialConditions
  {
    ASPECT_REGISTER_COMPOSITIONAL_INITIAL_CONDITIONS(SubdCompo,
                                                     "3D subduction",
                                                     "The compositions are based on the thickness of the oceanic "
                                                     "overriding, subducting and adjacent plates as specified by the plate "
                                                     "cooling model for the initial temperature. ")
  }
}
