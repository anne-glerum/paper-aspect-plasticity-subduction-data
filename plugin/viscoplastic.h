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


#ifndef __aspect__model_viscoplastic_h
#define __aspect__model_viscoplastic_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model in which all parameters depend on composition.
     * Furthermore, density is temperature dependent, and viscosity
     * depends on temperature, pressure and strain rate.
     *
     * Viscosity is computed based on diffusion and dislocation creep
     * and a Drucker-Prager yield criterion.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Viscoplastic : public MaterialModel::InterfaceCompatibility<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        /** This function computes which flow law results in the lowest viscosity,
         * indicating which deformation mechanism is dominant.
         * The meaning of the returned values:
         * -1 = diffusion creep
         * 0  = dislocation creep
         * 1  = plasticity
         */
        virtual double viscosity_ratio (const double temperature,
                                        const double pressure,
                                        const std::vector<double>    &compositional_fields,
                                        const SymmetricTensor<2,dim> &strain_rate,
                                        const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
        * Return true if the viscosity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the density() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the compressibility() function returns something that
        * may depend on the variable identifies by the argument.
        *
        * This function must return false for all possible arguments if the
        * is_compressible() function returns false.
        */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the specific_heat() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the thermal_conductivity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;
        /**
         * @}
         */

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
        /**
         * @}
         */


      private:

        /**
         * Function to calculate the compositional field with the
        *  highest value for infinity norm averaging of the viscosity.
         */
        int maximum_composition(const std::vector<double> &comp) const;

        /**
         * Function to calculate the harmonic average of the contributions
         * of each compositional field to the viscosity.
         */
        double harmonic_average(const std::vector<double> &comp,
                                const std::vector<double> &eta) const;

        /**
         * Function to calculate the arithmetic average of the contributions
         * of the compositional fields to the viscosity.
         */
        double arithmetic_average(const std::vector<double> &comp,
                                  const std::vector<double> &eta) const;

        /**
         * Function to calculate the geometric average of the contributions
         * of the compositional fields to the viscosity.
         */
        double geometric_average(const std::vector<double> &comp,
                                 const std::vector<double> &eta) const;

        /*
         * Function to calculate diffusion creep viscosity.
         */
        double diffusion(const double prefactor,
                         const double activation_energy,
                         const double activation_volume,
                         const double temperature,
                         const double pressure,
                         const double nu) const;

        /*
         * Function to calculate dislocation creep viscosity.
         */
        double dislocation(const double prefactor,
                           const double stress_exponent,
                           const double activation_energy,
                           const double activation_volume,
                           const double temperature,
                           const double pressure,
                           const double strain_rate,
                           const double nu) const;

        /*
         * Function to calculate plastic viscosity.
         */
        double plastic(const double phi,
                       const double cohesion,
                       const double pressure,
                       const double strain_rate) const;

        /*
         * Use a global initial strain rate instead of an initial
         * viscosity for each compositional field.
         */
        bool use_initial_strain_rate;

        /*
         * The value used as global initial strain rate.
         */
        double initial_strain_rate;

        /*
         * The number of compositional fields.
         */
        unsigned int                   n_compositional_fields;

        /*
         * List of angles of internal friction of the compositional fields.
         */
        std::vector<double>            phis_fields;

        /*
         * List of thermal conductivities of the compositional fields.
         */
        std::vector<double>            conductivities_fields;

        /*
         * List of heat capacities of the compositional fields.
         */
        std::vector<double>            capacities_fields;

        /*
         * List of reference densities of the compositional fields used in the density calculation.
         */
        std::vector<double>            refdens_fields;

        /*
         * List of reference temperatures of the compositional fields used in the density calculation.
         */
        std::vector<double>            reftemps_fields;

        /*
         * List of cohesions of the compositional fields used in the plastic viscosity calculation.
         */
        std::vector<double>            cohesions_fields;

        /*
         * List of strain rate exponents of the compositional fields used in the calculation of dislocation creep.
         */
        std::vector<double>            stress_exponents_fields;

        /*
         * List of activation volumes of the compositional fields used in the calculation of diffusion creep.
         */
        std::vector<double>            activation_volumes_diffusion_fields;

        /*
         * List of activation energies of the compositional fields used in the calculation of diffusion creep.
         */
        std::vector<double>            activation_energies_diffusion_fields;

        /*
         * List of prefactors of the compositional fields used in the calculation of diffusion creep.
         * This includes water fugacity, water fugaticty exponent, grain size and grain size exponent.
         */
        std::vector<double>            prefactors_diffusion_fields;

        /*
         * List of activation volumes of the compositional fields used in the calculation of dislocation creep.
         */
        std::vector<double>            activation_volumes_dislocation_fields;

        /*
         * List of activation energies of the compositional fields used in the calculation of dislocation creep.
         */
        std::vector<double>            activation_energies_dislocation_fields;

        /*
         * List of prefactors of the compositional fields used in the calculation of dislocation creep.
         * This includes water fugacity, water fugacity exponent, grain size and grain size exponent.
         */
        std::vector<double>            prefactors_dislocation_fields;

        /*
         * A scaling factor for dislocation creep viscosity.
         */
        std::vector<double>            nu_dislocation_fields;

        /*
         * A scaling factor for diffusion creep viscosity.
         */
        std::vector<double>            nu_diffusion_fields;

        /*
         * List of initial viscosities of the compositional fields
         * that are prescribed in the very first nonlinear iteration,
         * when not using an initial strain rate.
         */
        std::vector<double>            init_eta_fields;

        /**
         * Whether or not to correct creep prefactors for uniaxiality.
         * See for example the book of Gerya (2010), Section 6.2.
         */
        bool correct_uniaxiality;

        /**
         * The viscosity used for scaling the governing equations.
         */
        double reference_eta;

        /**
         * The viscosity prescribed during the very first nonlinear iteration in absence of compositional fields.
         */
        double initial_eta;

        /**
         * The minimum value of the viscosity cut-off used to restrict the overall viscosity ratio.
         */
        double minimum_eta;

        /**
         * The maximum value of the viscosity cut-off used to restrict the overall viscosity ratio.
         */
        double maximum_eta;

        /**
         * The type of averaging used for the contributions of the compositional fields to viscosity (and density).
         */
        std::string viscosity_averaging;

        /**
         * Whether or not to use strain rate weakening in the calculation of plastic viscosity.
         * If true, the plastic viscosity is multiplied with:
         * (std::max(1.0-(effective_strain_rate/reference_strain_rate),0.1))
         */
        bool strain_rate_weakening;

        /**
         * The reference strain rate used for strain rate weakening.
         */
        double ref_strain_rate;

        /**
         * Whether or not to use a harmonic average of the viscous and plastic viscosity
         * or to take the minimum of the two.
         */
        bool harmonic_plastic_viscous;

        /**
         * Whether or not to use a harmonic average of the effective and maximum viscosity
         * or to take the minimum. If true, the minimum viscosity is also added to the viscosity.
         */
        bool harmonic_max;

        /**
         * Whether or not the apply a function-specified weak zone as a last step in the viscosity calculation.
         */
        bool   weak_zone;

        /*
         * Parsed Function that specifies where to apply additional weakness
         * when weak_zone evaluates to true.
         * The weakened area is defined in Cartesian coordinates,
         * the weakening itself as a factor by which the viscosity is divided.
         */
        Functions::ParsedFunction<dim> weak_zone_function;

        /**
         * The constant thermal expansivity.
         */
        double thermal_alpha;

        /**
         * The thermal conductivity.
         */
        double k_value;

        /**
         * The reference density used for scaling the governing equations.
         * It is also used when no compositional fields are specified to compute
         * the density.
         */
        double reference_rho;

        /**
         * The reference specific heat used in calculating the reference thermal diffusivity.
         * It is also used when no compositional fields are specified to compute
         * the specific heat.
         */
        double reference_specific_heat;

        /**
         * The reference temperature used in the density calculation when no
         * compositional fields are specified.
         */
        double reference_T;

        /**
         * The activation energy used in the calculation of diffusion creep
         * when no compositional fields are specified.
         */
        double activation_energy_diffusion;

        /**
         * The activation volume used in the calculation of diffusion creep
         * when no compositional fields are specified.
         */
        double activation_volume_diffusion;

        /**
         * The prefactor used in the calculation of diffusion creep, including water fugacity,
         * water fugacity exponent, grain size and grain size exponent,
         * when no compositional fields are specified.
         */
        double prefactor_diffusion;

        /**
         * The activation energy used in the calculation of dislocation creep
         * when no compositional fields are specified.
         */
        double activation_energy_dislocation;

        /**
         * The activation volume used in the calculation of diffusion creep
         * when no compositional fields are specified.
         */
        double activation_volume_dislocation;

        /**
         * The prefactor used in the calculation of diffusion creep, including water fugacity,
         * water fugacity exponent, grain size and grain size exponent,
         * when no compositional fields are specified.
         */
        double prefactor_dislocation;

        /**
         * A scaling factor that can be used in the calculation of dislocation creep
         * when no compositional fields are specified.
         */
        double nu_dislocation;

        /**
         * A scaling factor that can be used in the calculation of diffusion creep
         * when no compositional fields are specified.
         */
        double nu_diffusion;

        /**
         * The exponent of strain rate used in the calculation of dislocation creep
         * when no compositional fields are specified. (For diffusion creep the
         * exponent is assumed to be 1.)
         */
        double stress_exponent;

        /**
         * The angle of internal friction used in the calculation of plasticity
         * when no compositional fields are specified.
         */
        double phi;

        /**
         * The cohesion used in the calculation of plasticity
         * when no compositional fields are specified.
         */
        double C;

    };

  }
}

#endif
