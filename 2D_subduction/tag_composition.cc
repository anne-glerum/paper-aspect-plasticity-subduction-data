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


#include "tag_composition.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    TagComposition<dim>::tag_additional_cells() const
    {
      if (this->get_dof_handler().n_dofs() != 0)
        {

          AssertThrow (this->n_compositional_fields() >= 1,
                       ExcMessage ("This refinement criterion can not be used when no "
                                   "compositional fields are active!"));

          QGauss<dim> quadrature (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);

          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature,
                                   update_quadrature_points | update_values);

          // the values of the compositional fields are stored as blockvectors for each field
          // we have to extract them in this structure
          std::vector<std::vector<double> > prelim_composition_values
          (this->n_compositional_fields(),
           std::vector<double> (quadrature.size()));

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();
          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                bool coarsen = false;
                bool refine = false;
                bool clear_refine = false;
                bool clear_coarsen = false;
                bool crust_present = false;
                bool slab_mantle_present = false;
                bool overriding_present = false;
                bool in_center_of_compo = false;

                fe_values.reinit(cell);

                for (unsigned int c = 0; c<this->n_compositional_fields(); c++)
                  {
                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                        prelim_composition_values[c]);
                  }


                for (unsigned int p=0; p<quadrature.size(); ++p)
                  {
                    if (prelim_composition_values[crust_refinement[0]][p] > 0.01)
                      {
                        crust_present = true;
                        //Crust will have smallest res, so not interested in other fields
                        break;
                      }
                    if (prelim_composition_values[slab_mantle_refinement[0]][p] > 0.1)
                      {
                        slab_mantle_present = true;
                        if (prelim_composition_values[slab_mantle_refinement[0]][p] >= 1.0)
                        {
                          in_center_of_compo = true;
                        }
                      }
                    if (prelim_composition_values[overriding_refinement[0]][p] > 0.1)
                      {
                        overriding_present = true;
                        if (prelim_composition_values[overriding_refinement[0]][p] >= 1.0)
                        {
                          in_center_of_compo = true;
                        }
                      }
                  }


                //Only continue if at least one is true

                    int maximum_refinement_level = max_level;
                    int minimum_refinement_level = min_level;

                    if (crust_present)
                      {
                        minimum_refinement_level = crust_refinement[1];
                        maximum_refinement_level = crust_refinement[2];
                      }
                    else if (slab_mantle_present)
                      {
                        minimum_refinement_level = slab_mantle_refinement[1];
                        maximum_refinement_level = slab_mantle_refinement[2];
                      }
                    else if (overriding_present)
                      {
                        minimum_refinement_level = overriding_refinement[1];
                        maximum_refinement_level = overriding_refinement[2];
                      }
                    else
                      {
                        minimum_refinement_level = mantle_refinement[0];
                        maximum_refinement_level = mantle_refinement[1];
                      }

                    const int cell_level = cell->level();
                    if (cell_level >= maximum_refinement_level)
                      {
                        clear_refine = true;
                      }
                    if (cell_level >  maximum_refinement_level)
                      {
                        coarsen = true;
                      }
                    if (cell_level <= minimum_refinement_level)
                      {
                        clear_coarsen = true;
                      }
                    if (cell_level < minimum_refinement_level)
                      {
                        refine = true;
                      }

                    if (clear_refine)
                      cell->clear_refine_flag ();
                    if (clear_coarsen)
                      cell->clear_coarsen_flag ();
                    if (refine)
                      cell->set_refine_flag ();
                    if (coarsen)
                      cell->set_coarsen_flag ();

              }
        }
    }

    template <int dim>
    void
    TagComposition<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition");
        {
          prm.declare_entry("Crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Slab mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the slab mantle, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Overriding plate refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the overriding plate, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the mantle, its minimum refinement level and "
                            "its maximum refinement level.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    TagComposition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        min_level = prm.get_integer("Minimum refinement level");
        max_level = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
        prm.enter_subsection("Composition");
        {

          const std::vector<int> crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Crust refinement")));

          crust_refinement = std::vector<unsigned int> (crust.begin(),crust.end());

          AssertThrow (crust_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (crust_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the crust cannot be "
                                   "greater than the maximum level of the whole model. "));

          const std::vector<int> slab_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Slab mantle refinement")));

          slab_mantle_refinement = std::vector<unsigned int> (slab_mantle.begin(),
                                                              slab_mantle.end());

          AssertThrow (slab_mantle_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (slab_mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (slab_mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the slab mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (slab_mantle_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the slab mantle cannot be "
                                   "greater than the maximum level of the whole model. "));


          const std::vector<int> overriding
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Overriding plate refinement")));

          overriding_refinement = std::vector <unsigned int> (overriding.begin(),
                                                              overriding.end());

          AssertThrow (overriding_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (overriding_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (overriding_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the overriding plate cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (overriding_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the overriding plate cannot be "
                                   "greater than the maximum level of the whole model. "));


          const std::vector<int> mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Mantle refinement")));

          mantle_refinement = std::vector <unsigned int> (mantle.begin(),
                                                          mantle.end());

          AssertThrow (mantle_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));
          
          AssertThrow (mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (mantle_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the mantle cannot be "
                                   "greater than the maximum level of the whole model. "));

          AssertThrow (crust_refinement[0] != slab_mantle_refinement[0] && \
                       crust_refinement[0] != overriding_refinement[0]  && \
                       crust_refinement[0] != mantle_refinement[0]  && \
                       slab_mantle_refinement[0] != overriding_refinement[0] && \
                       slab_mantle_refinement[0] != mantle_refinement[0] && \
                       overriding_refinement[0]  != mantle_refinement[0], 
                       ExcMessage ("Defined refinement fields the same. "));

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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(TagComposition,
                                              "tag composition",
                                              "A mesh refinement criterion that "
                                              "(de)flags cells for refinement and coarsening "
                                              "based on what composition is present. Different max "
                                              "and min refinement levels are set for the mantle, the crustal "
                                              "field, the slab mantle and the overriding plate. ")
  }
}
