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

#include "boxaddbi.h"

#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/base/std_cxx1x/bind.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    BoxAddBI<dim>::
    set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      int vertex_1, vertex_2, coord;
      // Horizontal domain boundaries
      if (BI == 2 || BI == 3)
        {
          vertex_1 = 2;
          vertex_2 = 3;
          coord    = 0;
        }

      // Vertical domain boundaries
      if (BI == 0 || BI == 1)
        {
          vertex_1 = 0;
          vertex_2 = 2;
          coord    = 1;
        }

      // Iterate over each cell to establish whether it lies on the specified boundary
      // and is part of the region that makes up the new boundary
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {
          if (cell->face(BI)->at_boundary())
            {
              if (cell->vertex(vertex_1)[coord] >= add_boundary_coord[0] &&
                  cell->vertex(vertex_2)[coord] <= add_boundary_coord[1])
                {
                  cell->face(BI)->set_boundary_id (2*dim);
                }
            }
        }
    }


    template <int dim>
    void
    BoxAddBI<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      AssertThrow(dim==2,ExcMessage("This geometry is currently only available for 2D."));

      GridGenerator::hyper_rectangle (coarse_grid,
                                      Point<dim>(),
                                      extents);

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        coarse_grid.begin_active()->face(f)->set_all_boundary_ids(f);

      //Make sure the right BI is set after refinement through the function set_boundary_indicators above
      coarse_grid.signals.post_refinement.connect
      (std_cxx1x::bind (&BoxAddBI<dim>::set_boundary_indicators,
                        std_cxx1x::cref(*this),
                        std_cxx1x::ref(coarse_grid)));

    }


    template <int dim>
    std::set<types::boundary_id>
    BoxAddBI<dim>::
    get_used_boundary_indicators () const
    {
      // Boundary indicators are zero through 2*dim
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim+1; ++i)
        {
          s.insert (i);
        }
      return s;
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    BoxAddBI<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("bottom", 2),
                  std::pair<std::string,types::boundary_id>("top",    3),
                  std::pair<std::string,types::boundary_id>("add_bi", 4)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("front",  2),
                  std::pair<std::string,types::boundary_id>("back",   3),
                  std::pair<std::string,types::boundary_id>("bottom", 4),
                  std::pair<std::string,types::boundary_id>("top",    5),
                  std::pair<std::string,types::boundary_id>("add_bi", 6)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
    BoxAddBI<dim>::
    get_periodic_boundary_pairs () const
    {
      // We don't need periodic boundaries,
      // so return empty set.
      std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundaries;
      return periodic_boundaries;
    }


    template <int dim>
    Point<dim>
    BoxAddBI<dim>::get_extents () const
    {
      return extents;
    }


    template <int dim>
    Point<dim>
    BoxAddBI<dim>::get_origin () const
    {
      // We don't need another origin than 0,0.
      Point<dim> box_origin;
      return box_origin;
    }


    template <int dim>
    double
    BoxAddBI<dim>::
    length_scale () const
    {
      return length_scale_factor*extents[0];
    }


    template <int dim>
    double
    BoxAddBI<dim>::depth(const Point<dim> &position) const
    {
      const double d = maximal_depth()-position(dim-1);

      return std::min (std::max (d, 0.), maximal_depth());
    }


    template <int dim>
    Point<dim>
    BoxAddBI<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point on the center axis of the domain
      Point<dim> p = extents/2;
      p[dim-1] = maximal_depth() - depth;
      return p;
    }


    template <int dim>
    double
    BoxAddBI<dim>::maximal_depth() const
    {
      return extents[dim-1];
    }

    template <int dim>
    bool
    BoxAddBI<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(this->get_free_surface_boundary_indicators().size() == 0 ||
                  this->get_timestep_number() == 0,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != 0,
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      for (unsigned int d = 0; d < dim; d++)
        if (point[d] > extents[d]+std::numeric_limits<double>::epsilon()*extents[d] ||
            point[d] < 0.-std::numeric_limits<double>::epsilon()*extents[d])
          return false;

      return true;
    }


    template <int dim>
    bool
    BoxAddBI<dim>::has_curved_elements() const
    {
      return false;
    }


    template <int dim>
    void
    BoxAddBI<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("BoxAddBI");
        {
          prm.declare_entry ("X extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in x-direction. Units: m.");
          prm.declare_entry ("Y extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in y-direction. Units: m.");
          prm.declare_entry ("Length scaling factor", "0.01",
                             Patterns::Double (0),
                             "A scaling factor with which the x-extent of the box will be "
                             "multiplied to obtain a length scale. Units: m.");
          prm.declare_entry ("Boundary indicator", "3",
                             Patterns::Integer (0),
                             "The boundary indicator number of the boundary of the box "
                             "geometry on which to specify an additional boundary indicator.");
          prm.declare_entry ("Additional boundary start", "0.4375",
                             Patterns::Double (0),
                             "Starting point of the X|Y range of the additional boundary. "
                             "Units: m.");
          prm.declare_entry ("Additional boundary end", "0.5625",
                             Patterns::Double (0),
                             "Starting point of the X|Y range of the additional boundary. "
                             "Units: m.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    BoxAddBI<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("BoxAddBI");
        {
          extents[0] = prm.get_double ("X extent");

          extents[1] = prm.get_double ("Y extent");

          length_scale_factor = prm.get_double ("Length scaling factor");

          BI = prm.get_integer ("Boundary indicator");

          add_boundary_coord[0] = prm.get_double ("Additional boundary start");
          add_boundary_coord[1] = prm.get_double ("Additional boundary end");
          Assert (add_boundary_coord[0] < add_boundary_coord[1],
                  ExcMessage ("Boundary start must be smaller than boundary end."));
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
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(BoxAddBI,
                                   "boxaddbi",
                                   "A 2d box geometry parallel to the coordinate directions. "
                                   "The extent of the box in each coordinate direction "
                                   "is set in the parameter file. The boxaddbi geometry labels its "
                                   "2*dim+1 sides as follows: boundary indicators 0 through 3 "
                                   "denote the left, right, bottom and top boundaries. "
                                   "On top of this is the boundary indicator 2*dim=4, "
                                   "set to a user-specified part of one of the boundaries. "
                                   " See also the documentation of the deal.II class "
                                   "``GeometryInfo''.")
  }
}
