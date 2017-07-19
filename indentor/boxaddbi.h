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



#ifndef __aspect__geometry_model_boxaddbi_h
#define __aspect__geometry_model_boxaddbi_h

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a 2D box geometry of certain width and height.
     * An additional boundary indicator can be specified, such that different
     * boundary conditions can be applied to one boundary of the domain.
     * In particular, this is useful for the punch benchmark.
     */
    template <int dim>
    class BoxAddBI : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Return a point that denotes the size of the box in each dimension
         * of the domain.
         */
        virtual
        Point<dim> get_extents () const;

        /**
         * Return a point that denotes the lower left corner of the box
         * domain.
         */
        virtual
        Point<dim> get_origin () const;

        /**
         * Return the typical length scale one would expect of features in this geometry,
         * assuming realistic parameters.
         *
         * We return 1/100th of the diameter of the box.
         */
        virtual
        double length_scale () const;

        /**
         * Return the depth that corresponds to the given
         * position. The documentation of the base class (see
         * GeometryModel::Interface::depth()) describes in detail how
         * "depth" is interpreted in general.
         *
         * Computing a depth requires a geometry model to define a
         * "vertical" direction. The current class considers the
         * $(0,1)^T$ vector in 2d (and the $(0,0,1)^T$ vector in 3d)
         * as vertical and considers the "top" boundary as the
         * "surface". In almost all cases one will use a gravity model
         * that also matches these definitions.
         */
        virtual
        double depth(const Point<dim> &position) const;

        virtual
        Point<dim> representative_point(const double depth) const;

        virtual
        double maximal_depth() const;

        /**
         * Return the set of boundary indicators that are used by this model. This
         * information is used to determine what boundary indicators can be used in
         * the input file.
         *
         * The boxaddbi model uses boundary indicators zero through 2*dim=4, with the first
         * two being the faces perpendicular to the x-axis, the next two perpendicular
         * to the y-axis and 42 a userspecified subarea of one of the preceding boundaries.
         */
        virtual
        std::set<types::boundary_id>
        get_used_boundary_indicators () const;

        /**
         * Return a mapping from symbolic names of each part of the boundary
         * to the corresponding boundary indicator. This allows users to
         * specify *names*, not just *numbers* in their input files when
         * describing which parts of the boundary have to satisfy which
         * boundary conditions.
         *
         * This geometry returns the map <code>{{"left"->0}, {"right"->1},
         * {"bottom"->2}, {"top"->3}, {"add_bi"->4}}</code> in 2d, and <code>{{"left"->0},
         * {"right"->1}, {"front"->2}, {"back"->3}, {"bottom"->4},
         * {"top"->5},{"add_bi"->6}</code> in 3d.
         */
        virtual
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const;

        /**
         * Return the set of periodic boundaries as described in the input
         * file.
         */
        virtual
        std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
        get_periodic_boundary_pairs () const;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        virtual
        bool
        point_is_in_domain(const Point<dim> &p) const;

        /**
         * @copydoc Interface::has_curved_elements()
         *
         * A box has only straight boundaries and cells, so return false.
         */
        virtual
        bool
        has_curved_elements() const;

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
         * Bind boundary indicators to child cells after each mesh refinement round.
         */
        virtual
        void
        set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const;

        /**
         * Extent of the box in x- and y-direction.
         */
        Point<dim> extents;

        /**
         * The factor to scale the x-extent of the box with to obtain a length scale.
         */
        double length_scale_factor;

        /**
         * The boundary indicator of the existing domain boundary on which the new indicator should lie.
         */
        int BI;

        /**
         * The begin and end x or y value of the stretch of boundary with the new boundary indicator.
         */
        double add_boundary_coord[2];

    };
  }
}


#endif
