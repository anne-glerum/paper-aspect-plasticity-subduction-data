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


#include "2Dsubd_compo.h"
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
      // Box. Verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      Assert (geometry != 0,
              ExcMessage ("This initial condition can only be used if the geometry "
                          "is a box."));

      // Depth of the current point
      const double depth = geometry->depth(position);
      // Uniform thermal diffusivity
      const double thermal_diffusivity = 1.0e-6;
      // The height of the domain and of the maximum plate thickness
       const double domain_width = geometry->get_extents()[0];

      // Determine plate ages based on distance from left or right domain boundary
      // Subducting plate (SP) starts from the right, overriding plate (OP) from the left
      // Leave 5 km free to ensure decoupling from side boundaries
      const double age_OP = (position[0]<5000.0) ? 0.0 : (position[0]-5000.0) * age_OP_max / x_trench;
      const double age_SP = (position[0]>domain_width-5000.0)? 0.0 : (domain_width-position[0]-5000.0) * age_SP_max / (domain_width - x_trench);

      // Horizontal (x) position of the slab tip for slab extending
      // 206 km into mantle
      const double slab_tip = x_trench - std::cos(slab_dip) * 206000.0;

      // The SP
      if (position[0] >= x_trench)
         { 
          // Crustal layer
    	  // Starts 21 km after mantle layer ends
          if (depth <= d_crust && \
              depth <= (2.32*std::sqrt(thermal_diffusivity*age_SP))-20800.0)
             return (n_comp == 3 ) ? 1.0 : 0.0;
          // Mantle layer
          // If we're at the vertical domain wall, don't set composition
          else if (depth <= (2.32*std::sqrt(thermal_diffusivity*age_SP)))
             return (n_comp == 2 ) ? (age_SP == 0.0 ? 0.0 : 1.0) : 0.0;
          else 
          return (n_comp == 0) ? 1.0 : 0.0;
         }
      // The slab tip
      else if (position[0] >= slab_tip && \
               position[0] < x_trench &&  \
               depth >= std::tan(slab_dip) * (x_trench-position[0]))
              {
               //Prescribe T in slab tip by transposing the top of the slab to the surface
               Point<dim> new_position; 
               new_position[0] = position[0];
               new_position[dim-1] = position[dim-1]+std::tan(slab_dip) * (x_trench-position[0]);
               const double new_depth = geometry->depth(new_position);

               if (new_depth <= d_crust)
                  return (n_comp == 3 ) ? 2.0 : 0.0;
               else if (new_depth > d_crust && \
                        new_depth <= (2.32*std::sqrt(thermal_diffusivity*age_SP)) )
                  return (n_comp == 2) ? 1.0 : 0.0;
               else
                  return (n_comp == 0) ? 1.0 : 0.0;
              }
      // The OP
      else if (depth < std::tan(slab_dip) * (x_trench-position[0]) && \
               depth <= (2.32*std::sqrt(thermal_diffusivity*age_OP)))
               {
                  return (n_comp == 1 ) ? (age_OP == 0.0 ? 0.0 : 1.0) : 0.0;
              }
      else
       return (n_comp == 0) ? 1.0 : 0.0;


    }

    template <int dim>
    void
    SubdCompo<dim>::declare_parameters (ParameterHandler &)
    {
    }


    template <int dim>
    void
    SubdCompo<dim>::parse_parameters (ParameterHandler &prm)
    {
      const bool use_years_instead_of_seconds = prm.get_bool ("Use years in output instead of seconds");
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plate cooling");
        {
          age_OP_max = (use_years_instead_of_seconds ? year_in_seconds : 1.0 ) * prm.get_double ("Overriding plate age at trench");
          age_SP_max = (use_years_instead_of_seconds ? year_in_seconds : 1.0 ) * prm.get_double ("Subducting plate age at trench");
          x_trench = prm.get_double ("Horizontal trench position");
          d_max = prm.get_double ("Maximum oceanic plate thickness");
          d_crust = prm.get_double ("Crustal thickness");
          slab_dip = prm.get_double ("Slab dip angle")*numbers::PI/180.0;
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
                                                     "2D subduction",
                                                     "The compositions are based on the thickness of the oceanic "
                                                     "overriding and subducting plates as specified by the plate "
                                                     "cooling model initial temperature conditions. ")
  }
}
