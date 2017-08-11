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


#include "3Dsubd_temp.h"
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


      //Get the (adiabatic) temperature at the top and bottom boundary of the model
      const double Ts = T_boundary->minimal_temperature(this->get_fixed_temperature_boundary_indicators());
      const double Tas = this->get_adiabatic_surface_temperature();
      const double Tb = T_boundary->maximal_temperature(this->get_fixed_temperature_boundary_indicators());

      const double depth = geometry->depth(position);
      double kap = 1.0e-6;

      //Look up material properties for calculation of thermal diffusivity
      //Determine plate age based on distance from domain boundary + 5 km
      //const double age_OP = position[0]-5000.0 * (1.0/v_spread_OP);
      //const double age_SP = (geometry->get_extents()[0]-position[0]+5000.0) * (1.0/v_spread_SP);
      const double age_OP = (position[0]<5000.0) ? 0.0 : (position[0]-5000.0) * (1.0/v_spread_OP);
      const double age_SP = (position[0]>geometry->get_extents()[0]-5000.0)? 0.0 : (geometry->get_extents()[0]-position[0]-5000.0) * (1.0/v_spread_SP);

      // The height of the domain and of the maximum plate thickness
      const double domain_height = geometry->get_extents()[dim-1];
      const double z_max = domain_height-d_max;

      // The slab dip angle in radians
      const double slab_angle = 29.0 * numbers::PI / 180.0;

      // Trench position
      const double trench = 1600000.0;

      // Horizontal (x) position of the slab tip
      const double slab_tip = trench - std::cos(slab_angle) * 206000.0;

      // The parameters needed for the plate cooling temperature calculation
      double temp = 0.0;
      const int n_sum = 100;
      double sum_OP = 0.0;
      double sum_SP = 0.0;
      double sum_AP = 0.0;

       
      for (int i=1;i<=n_sum;i++)
           {
            sum_OP += (1.0/i) * 
                          (exp((-kap*i*i*numbers::PI*numbers::PI*age_OP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*depth/d_max));

            sum_SP += (1.0/i) * 
                          (exp((-kap*i*i*numbers::PI*numbers::PI*age_SP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*depth/d_max));

            sum_AP += (1.0/i) * 
                          (exp((-kap*i*i*numbers::PI*numbers::PI*age_AP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*depth/d_max));
           }

      // The SP
      if (position[0] >= trench && position[dim-1] >= z_max && position[dim-2] > plate_width) 
         { 
           temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_SP));
         }
      // The slab tip
      else if (position[0] >= slab_tip && \
               position[0] < trench &&  \
               position[dim-1] <= domain_height-std::tan(slab_angle) * (trench-position[0]) && \
               position[dim-1] >= z_max-std::tan(slab_angle) * (trench-position[0]) && \
               position[dim-2] > plate_width)
              {
               //Prescribe T in slab tip by transposing the top of the slab to the surface
               Point<dim> new_position; 
               new_position[0] = position[0];
               new_position[dim-1] = position[dim-1]+std::tan(slab_angle) * (trench-position[0]);
               const double new_depth = geometry->depth(new_position);

               sum_SP = 0.0;
               for (int i=1;i<=n_sum;i++)
               sum_SP += (1.0/i) * 
                          (exp((-kap*i*i*numbers::PI*numbers::PI*age_SP)/(d_max*d_max)))*
                          (sin(i*numbers::PI*new_depth/d_max));

               temp = std::min(Tm,Ts + (Tm - Ts) * ((new_depth / d_max) + (2.0 / numbers::PI) * sum_SP));
              }
      // The OP
      else if (position[dim-1] >= z_max && \
               position[dim-1] > domain_height-std::tan(slab_angle) * (trench-position[0]) && \
               position[dim-2] > plate_width)
              {
               temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_OP));
              }
      // The AP
      else if (position[dim-2] <= plate_width - 40000.0 && position[dim-1] >= z_max)
      {
         temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_AP));
      }
      // The WZ
      else if (position[dim-2] > plate_width - 40000.0 && position[dim-2] <= plate_width && position[dim-1] >= z_max)
      {
       if (position[0] >= trench )
         temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_SP)); 
       else 
         temp = std::min(Tm,Ts + (Tm - Ts) * ((depth / d_max) + (2.0 / numbers::PI) * sum_OP));
      }
      // The mantle: linear gradient
      else
          {
           temp = -((Tb - Tm) / ((domain_height - d_max) / 1000.0)) * (position[dim-1] / 1000.0) + Tb; 
          }

      if(temp <=0.0 || temp >= 2000.0)
      std::cout << "T " << temp << " " << position << std::endl;
      return std::max(0.0,std::min(temp,2000.0));
    }

    template <int dim>
    void
    SubdTemp<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plate cooling");
        {
          prm.declare_entry ("Spreading velocity subducting plate", "0.0407653503",
                             Patterns::Double (0),
                             "The spreading velocity of the subducting plate, used for the calculation "
                             "of the plate cooling model temperature. This represents the right plate. Units: m/years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "m/seconds otherwise.");
          prm.declare_entry ("Spreading velocity overriding plate", "0.1087076011",
                             Patterns::Double (0),
                             "The spreading velocity of the overriding plate, used for the calculation "
                             "of the plate cooling model temperature. This represents the right plate. Units: m/years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "m/seconds otherwise.");
          prm.declare_entry ("Width of subducting and overriding plate", "400000.0",
                             Patterns::Double (0),
                             "The width of the subducting and overriding plate. "
                             "The rest of the domain in the y-direction is adjacent plate. "
                             "Units: m. ");
          prm.declare_entry ("Maximum oceanic plate thickness", "125000.0",
                             Patterns::Double (0),
                             "The maximum thickness of an oceanic plate in the plate cooling model "
                             "for when time goes to infinity. Units: m. " );
          prm.declare_entry ("Maximum oceanic plate temperature", "1593.00",
                             Patterns::Double (0),
                             "The maximum temperature of an oceanic plate in the plate cooling model "
                             "for when time goes to infinity. Units: K. " );
          prm.declare_entry ("Age adjacent plate", "58.87352819e6",
                             Patterns::Double (0),
                             "The age of the adjacent plate that is fixed with distance "
                             "Units: years or seconds. " );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
 

    template <int dim>
    void
    SubdTemp<dim>::parse_parameters (ParameterHandler &prm)
    {
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
          Tm = prm.get_double ("Maximum oceanic plate temperature");
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
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(SubdTemp,
                                       "3D subduction",
                                       "An initial temperature field determined from the plate"
                                       "cooling model for a 3D subduction setup with a subducting "
                                       "plate, overriding plate, adjacent plate and the inbetween transform zone. ")
  }
}
