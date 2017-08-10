/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
/*  $Id$  */


#include "shear_band_angle.h"
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/mpi.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ShearBandAngle<dim>::execute (TableHandler &statistics)
    {

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities)
                                            .degree+1);

      const unsigned int n_q_points = quadrature_formula.size();

      std::vector<Point<dim> > position_values(n_q_points);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      // xposition: where measurements should be taken
      // yposition: the height of the maximum strain rate at xposition
      // max_strain_rate: the maximum strain rate found at xposition
      const double xposition[4] = {17400.0,19400.0,20600.0,22600.0};
      double yposition[4] = {0.0,0.0,0.0,0.0};
      double max_strain_rate[4] = {0.0,0.0,0.0,0.0};

      // Use this as a container for the strain rate
      typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                     this->n_compositional_fields());

      // Iterate over all cells to check whether they are at the xposition and then find
      // the maximum strain rate
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // retrieve quadrature points locations
            position_values = fe_values.get_quadrature_points();

            // retrieve the strain rate
            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                in.strain_rate);

            // determine bandwidth in which to search for the maximum strain rate
            const double delta_x = 0.5 * (position_values[1][0]-position_values[0][0]);

            // find the maximum strain rate at each xposition
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
            	for (unsigned int p = 0; p < 4; ++p)
            	{
            		if (position_values[q][0] >= (xposition[p]-delta_x) &&
            		    position_values[q][0] <= (xposition[p]+delta_x) &&
            		    in.strain_rate[q].norm() > max_strain_rate[p])
            	    {
                	    max_strain_rate[p] = in.strain_rate[q].norm();
            		    yposition[p] = position_values[q][1];
                    }
            	}
              }
          }

         // a struct for maximum strain rate and rank of the current processor
         struct {double strainrate_max; int rank;} local[4], global[4];

         // obtain rank of current processor
         int myrank;
         const int root = 0;
         MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

         // fill struct with max strain rate and rank of current processor
         for (unsigned int p = 0; p < 4; ++p)
     	 {
        	 local[p].strainrate_max = std::log10(std::max(1e-20,max_strain_rate[p]));
        	 local[p].rank = myrank;
     	 }

         // find the processor with the highest maximum strain rate at xposition
         MPI_Allreduce(&local,&global,4,MPI_DOUBLE_INT,MPI_MAXLOC, this->get_mpi_communicator());

         // collect data of processor with highest maximum strain rate on root
         double global_yposition[4] = {0.0, 0.0, 0.0, 0.0};

         for (unsigned int p = 0; p < 4; ++p)
         {
        	 if (myrank == root)
        	 {
        		 MPI_Status status[4];
        		 if (myrank == global[p].rank)
        		 {
        			 global_yposition[p] = yposition[p];
        		 }
        		 else
        		 {
        			 const unsigned int mpi_tag = 110 * (p+1);
        			 MPI_Recv (&global_yposition[p],1, MPI_DOUBLE, global[p].rank, mpi_tag, this->get_mpi_communicator(), &status[p]);
        		 }
        	 }
        	 if (myrank != root && myrank == global[p].rank)
        	 {
        		 const unsigned int mpi_tag = 110 * (p+1);
        		 MPI_Send(&yposition[p], 1, MPI_DOUBLE, root, mpi_tag, this->get_mpi_communicator());
        	 }

         }

         if (myrank == root)
         {
         // calculate the average shear band angle on root
         const double dx = 2000.0;
         const double dy_left = global_yposition[0]-global_yposition[1];
         const double dy_right = global_yposition[3] - global_yposition[2];
         const double right_angle = atan(dy_right/dx)*180.0/numbers::PI;
         const double left_angle = atan(dy_left/dx)*180.0/numbers::PI;
         const double shear_band_angle = (right_angle + left_angle)/2.0;

          // fill statistics file
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.add_value ("Average shear band angle (degrees)", shear_band_angle);
          statistics.set_precision ("Average shear band angle (degrees)", 8);
          statistics.set_scientific ("Average shear band angle (degrees)", true);


      std::ostringstream output;
      output.precision(3);

        output << shear_band_angle
               << " degrees";

      return std::pair<std::string, std::string> ("Average shear band angle:",
                                                  output.str());
         }
      return std::pair<std::string, std::string> ("Calculating average shear band angle:",
              "please wait");
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ShearBandAngle,
                                  "shear band angle",
                                  "A postprocessor that computes the average shear band angle "
                                  "for the Kaus 2010 brick benchmark.  "
                                  "At 4 given locations we look for the maximum strain rate and the "
                                  "corresponding height to calculate the angle. "
                                 )
  }
}
