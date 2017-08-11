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


#include "2Dsubd_temp.h"
#include <aspect/geometry_model/box.h>
#include <aspect/boundary_temperature/box.h>


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    SubdTemp<dim>::
    initial_temperature (const Point<dim> &position) const
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
    
      //Get the (adiabatic) temperature at the top and bottom boundary of the model
      const double Ts = T_boundary->minimal_temperature(this->get_fixed_temperature_boundary_indicators());
      const double Tas = this->get_adiabatic_surface_temperature();
      const double Tb = T_boundary->maximal_temperature(this->get_fixed_temperature_boundary_indicators());

      const double depth = geometry->depth(position);
      const double thermal_diffusivity = 1.0e-6;
      // The height of the domain and of the maximum plate thickness
      const double domain_height = geometry->get_extents()[dim-1];
      const double z_max = domain_height-d_max;
      const double domain_width = geometry->get_extents()[0];

      //Determine plate age based on distance from domain boundary + 5 km
      //const double age_OP = position[0]-5000.0 * (1.0/v_spread_OP);
      //const double age_SP = (geometry->get_extents()[0]-position[0]+5000.0) * (1.0/v_spread_SP);
      const double age_OP = (position[0]<5000.0) ? 0.0 : (position[0]-5000.0) * age_OP_max / x_trench;
      const double age_SP = (position[0]>domain_width-5000.0)? 0.0 : (domain_width-position[0]-5000.0) * age_SP_max / (domain_width - x_trench);

      // Horizontal (x) position of the slab tip
      const double slab_tip = x_trench - std::cos(slab_dip) * 206000.0;

      // The parameters needed for the plate cooling temperature calculation
      double temp = 0.0;
      const int n_sum = 100;
      double sum_OP = 0.0;
      double sum_SP = 0.0;

      for (int i=1;i<=n_sum;i++)
           {
            sum_OP += (1.0/i) * 
                          (exp((-thermal_diffusivity*i*i*numbers::PI*numbers::PI*age_OP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*depth/d_max));

            sum_SP += (1.0/i) * 
                          (exp((-thermal_diffusivity*i*i*numbers::PI*numbers::PI*age_SP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*depth/d_max));
           }

      // The SP
      if (position[0] >= x_trench && position[dim-1] >= z_max)
         { 
           temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_SP));
         }
      // The slab tip
      else if (position[0] >= slab_tip && \
               position[0] < x_trench &&  \
               depth >= std::tan(slab_dip) * (x_trench-position[0]) && \
               depth <= domain_height - z_max + std::tan(slab_dip) * (x_trench-position[0]))
              {
               //Prescribe T in slab tip by transposing the top of the slab to the surface
               Point<dim> new_position; 
               new_position[0] = position[0];
               new_position[dim-1] = position[dim-1]+std::tan(slab_dip) * (x_trench-position[0]);
               const double new_depth = geometry->depth(new_position);

               sum_SP = 0.0;
               for (int i=1;i<=n_sum;i++)
               sum_SP += (1.0/i) * 
                          (exp((-thermal_diffusivity*i*i*numbers::PI*numbers::PI*age_SP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*new_depth/d_max));

               temp = std::min(Tm,Ts + (Tm - Ts) * ((new_depth / d_max) + (2.0 / numbers::PI) * sum_SP));
              }
      // The OP
      else if (position[dim-1] >= z_max && \
               depth < std::tan(slab_dip) * (x_trench-position[0]))
              {
               temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_OP));
              }
      // The mantle: linear gradient
      else
          {
           temp = -((Tb - Tm) / (z_max / 1000.0)) * (position[dim-1] / 1000.0) + Tb;
          }

      return temp;
    }

    template <int dim>
    void
    SubdTemp<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plate cooling");
        {
          prm.declare_entry ("Subducting plate age at trench", "100e6",
                             Patterns::Double (0),
                             "The age of the subducting plate, used for the calculation "
                             "of the plate cooling model temperature. This represents the right plate. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Overriding plate age at trench", "20e6",
                             Patterns::Double (0),
                             "The spreading velocity of the overriding plate, used for the calculation "
                             "of the plate cooling model temperature. This represents the left plate. Units: m/years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "m/seconds otherwise.");
          prm.declare_entry ("Maximum oceanic plate thickness", "125000.0",
                             Patterns::Double (0),
                             "The maximum thickness of an oceanic plate in the plate cooling model "
                             "for when time goes to infinity. Units: m. " );
          prm.declare_entry ("Maximum oceanic plate temperature", "1593.00",
                             Patterns::Double (0),
                             "The maximum temperature of an oceanic plate in the plate cooling model "
                             "for when time goes to infinity. Units: K. " );
          prm.declare_entry ("Crustal thickness", "10400.0",
                             Patterns::Double (0),
                             "The thickness of the crust of the subducting plate. "
                             "Units: m. " );
          prm.declare_entry ("Horizontal trench position", "1600000.0",
                             Patterns::Double (0),
                             "The horizontal position of the trench at the surface. Here the overriding plate ends. "
                             "Units: m. " );
          prm.declare_entry ("Slab dip angle", "29.0",
                             Patterns::Double (0),
                             "The uniform dip of the subduction slab. "
                             "Units: degrees. " );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
 

    template <int dim>
    void
    SubdTemp<dim>::parse_parameters (ParameterHandler &prm)
    {
      const bool use_years_instead_of_seconds = prm.get_bool ("Use years in output instead of seconds");
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plate cooling");
        {
          age_OP_max = (use_years_instead_of_seconds ? year_in_seconds : 1.0 ) * prm.get_double ("Overriding plate age at trench");
          age_SP_max = (use_years_instead_of_seconds ? year_in_seconds : 1.0 ) * prm.get_double ("Subducting plate age at trench");
          d_max = prm.get_double ("Maximum oceanic plate thickness");
          Tm = prm.get_double ("Maximum oceanic plate temperature");
          x_trench = prm.get_double ("Horizontal trench position");
          d_crust = prm.get_double ("Crustal thickness");
          slab_dip = prm.get_double ("Slab dip angle")*numbers::PI/180.0;
          AssertThrow (slab_dip < numbers::PI * 0.5, ExcMessage("Slab dip angle too great for plugin."));
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
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(SubdTemp,
                                       "2D subduction",
                                       "An initial temperature field determined from the plate"
                                       "cooling model, including a subducting and an overriding plate. ")
  }
}
