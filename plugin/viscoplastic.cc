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


#include "viscoplastic.h"
#include <deal.II/base/parameter_handler.h>


using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Get the field nr. of the field with the maximum value.
     */
    template <int dim>
    int
    Viscoplastic<dim>::
    maximum_composition(const std::vector<double> &comp) const
    {
      int maximum_comp = 0;
      double max = comp[0];

      for (unsigned int i = 1; i<comp.size(); i++)
        {
          if (comp[i] > max)
            {
              maximum_comp = i;
              max = comp[i];
            }
        }

      return maximum_comp;

    }

    /*
     * Compute the harmonic average of the compositional
     * contributions.
     */
    template <int dim>
    double
    Viscoplastic<dim>::
    harmonic_average(const std::vector<double> &comp,
                     const std::vector<double> &eta) const
    {
      double visc = 0.0;
      double compo_total = 0.0;

      for (unsigned int i = 0; i<n_compositional_fields; ++i)
        {
          visc += std::max(std::min(comp[i],1.0),0.0) / eta[i];
          compo_total += std::max(std::min(comp[i],1.0),0.0);
        }

      visc /= compo_total;
      visc = 1.0 / visc;

      return visc;
    }

    /*
     * Compute the arithmetic average of the compositional
     * contributions.
     */
    template <int dim>
    double
    Viscoplastic<dim>::
    arithmetic_average(const std::vector<double> &comp,
                       const std::vector<double> &eta) const
    {
      double visc = 0.0;
      double compo_total = 0.0;

      for (unsigned int n=0; n<n_compositional_fields; n++)
        {
          visc += std::max(std::min(comp[n],1.0),0.0) * eta[n];
          compo_total += std::max(std::min(comp[n],1.0),0.0);
        }

      visc /= compo_total;

      return visc;
    }

    /*
     * Compute the geometric average of the compositional
     * contributions.
     */
    template <int dim>
    double
    Viscoplastic<dim>::
    geometric_average(const std::vector<double> &comp,
                      const std::vector<double> &eta) const
    {
      double visc = 0.0;
      double compo_total = 0.0;

      for (unsigned int n=0; n<n_compositional_fields; n++)
        {
          visc += std::max(std::min(comp[n],1.0),0.0) * log10(eta[n]);
          compo_total += std::max(std::min(comp[n],1.0),0.0);
        }

      visc /= compo_total;
      visc = pow(10,visc);

      return visc;
    }


    template <int dim>
    double
    Viscoplastic<dim>::
    diffusion(const double prefactor,
              const double activation_energy,
              const double activation_volume,
              const double temperature,
              const double pressure,
              const double nu) const
    {
      // The gas constant
      const double R = 8.314;
      // Equation (6). of Glerum et al. (2017).
      // Assume n = 1 for diffusion creep.
      // Cut off temperature at 1 K to prevent division by zero
      // (for example in benchmark cases)
      // viscosity = 0.5 * beta / B * exp((Q+PV)/RT)
      return 0.5 * nu * (1e0/prefactor)*exp((activation_energy+activation_volume*pressure)/(R*std::max(1.0,temperature)));
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    dislocation(const double prefactor,
                const double stress_exponent,
                const double activation_energy,
                const double activation_volume,
                const double temperature,
                const double pressure,
                const double strain_rate_norm,
                const double nu) const
    {
      // The gas constant
      const double R = 8.314;
      // Equation (6). of Glerum et al. (2017).
      // Cut off temperature at 1 K to prevent division by zero
      // (for example in benchmark cases)
      // viscosity = 0.5 * beta * B**(-1/n) * effective_strainrate**((1-n)/n) * exp((Q+PV)/nRT)
      return (0.5 * nu * std::pow(prefactor,-1e0/stress_exponent)*
              std::pow(strain_rate_norm,(1e0-stress_exponent)/stress_exponent)*
              exp((activation_energy+activation_volume*pressure)/(stress_exponent*R*std::max(1.0,temperature))));
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    plastic(const double phi,
            const double cohesion,
            const double pressure,
            const double strain_rate_norm) const

    {
      // The yield strength
      double strength = 0;
      // In 3D, we approximate the hexagonal Mohr-Coulomb yield surface
      // with a smooth Drucker-Prager cone that matches the Mohr-Coulomb
      // surface either at the outer edges (ct = -1) or the inner edges
      // (ct = +1). Tensile stresses (negative) are cut off at zero.
      // Equation (9) of Glerum et al. (2017).
      if (dim==3)
        {
          const double ct = -1.; //1.;
          strength = ((6.0*cohesion*std::cos(phi))/(std::sqrt(3.0)*(3.0+ct*std::sin(phi)))) +
                     ((6.0*std::sin(phi))/(std::sqrt(3.0)*(3.0+ct*std::sin(phi)))) * std::max(pressure,0.0);
        }
      // In 2D incompressible plane strain,
      // we match the Mohr-Coulomb surface and the Drucker-Prager surface:
      // Equation (8) of Glerum et al. (2017).
      else strength = std::max(pressure,0.0) * std::sin(phi) + cohesion * std::cos(phi);
      return strength / (2.0*strain_rate_norm);
    }


    template <int dim>
    double
    Viscoplastic<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position) const
    {
      // The holder for plastic viscosity.
      double viscosity_plastic = 0.0;
      // The holder for viscous viscosity.
      double viscosity_viscous = 0.0;
      // The holder for the final viscosity.
      double viscosity = 0.0;

      // The very first nonlinear iteration (during the first time step,
      // when the initial guess for velocity/strain rate is zero),
      // we either use constant initial viscosities per material
      // or a global initial strain rate.
      if (this->get_timestep_number() == 0 && strain_rate.norm() == 0. && !use_initial_strain_rate)
        {
          // Multiple compositions
          if (n_compositional_fields>0)
            {
              if (viscosity_averaging=="Max")
                {
                  int max_comp = maximum_composition(composition);
                  viscosity = init_eta_fields[max_comp];
                }
              else if (viscosity_averaging=="Harmonic")
                {
                  viscosity = harmonic_average(composition, init_eta_fields);
                }
              else if (viscosity_averaging=="Arithmetic")
                {
                  viscosity = arithmetic_average(composition, init_eta_fields);
                }
              else
                {
                  viscosity = geometric_average(composition, init_eta_fields);
                }
            }
          // No compositional fields
          else
            viscosity = initial_eta;

          if (weak_zone)
            viscosity /= weak_zone_function.value(position);

          if (harmonic_max)
            {
              viscosity = 1.0 / ((1.0 / viscosity) + (1.0 / maximum_eta));
              viscosity += minimum_eta;
            }
          else
            viscosity = std::max(std::min(viscosity,maximum_eta),minimum_eta);

          return viscosity;
        }

      // Use the initial strain rate during the very first nonlinear iteration,
      // or compute the viscosity in the regular way.

      // Calculate the second invariant of the deviatoric strain rate tensor
      // The effective strain rate is sqrt(0.5*strain_rate_deviator*strain_rate_deviator).
      const SymmetricTensor<2,dim> strain_rate_dev = deviator(strain_rate);
      const double strain_rate_dev_inv = (use_initial_strain_rate && this->get_timestep_number() == 0 && strain_rate.norm() == 0.) ?
                                         initial_strain_rate :
                                         std::sqrt(0.5) * strain_rate_dev.norm();

      // Multiple compositional fields
      if (composition.size()>0)
        {
          if (viscosity_averaging=="Max")
            {
              // The field with the highest value
              int max_comp = maximum_composition(composition);
              // Compute the inverse of the diffusion viscosity.
              const double visc_diffusion_inverse = 1.0 / diffusion(prefactors_diffusion_fields[max_comp],
                                                                    activation_energies_diffusion_fields[max_comp],
                                                                    activation_volumes_diffusion_fields[max_comp],
                                                                    temperature,
                                                                    pressure,
                                                                    nu_diffusion_fields[max_comp]);

              // Compute the inverse of the dislocation viscosity.
              const double visc_dislocation_inverse = 1.0 / dislocation(prefactors_dislocation_fields[max_comp],
                                                                        stress_exponents_fields[max_comp],
                                                                        activation_energies_dislocation_fields[max_comp],
                                                                        activation_volumes_dislocation_fields[max_comp],
                                                                        temperature,
                                                                        pressure,
                                                                        strain_rate_dev_inv,
                                                                        nu_dislocation_fields[max_comp]);

              // Compute the composite viscosity.
              viscosity_viscous = 1.0 / (visc_diffusion_inverse + visc_dislocation_inverse);

              // Compute the plastic viscosity.
              viscosity_plastic =  plastic(phis_fields[max_comp],
                                           cohesions_fields[max_comp],
                                           pressure,
                                           strain_rate_dev_inv);

              // Apply strain rate weakening of plastic viscosity if needed.
              if (strain_rate_weakening)
                viscosity_plastic *= (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1));

              // Either average the plastic and viscous viscosity harmonically,
              // or take the minimum. Equation (12).
              if (harmonic_plastic_viscous)
                viscosity = 1.0 / ((1.0/viscosity_viscous) + (1.0/viscosity_plastic));
              else
                viscosity = std::min(viscosity_viscous,viscosity_plastic);
            }
          // Harmonic, arithmetic or geometric averaging
          else
            {
              // The viscosities for each compositional field.
              std::vector<double> visc_diffusion_inverse(n_compositional_fields);
              std::vector<double> visc_dislocation_inverse(n_compositional_fields);
              std::vector<double> visc_viscous(n_compositional_fields);
              std::vector<double> visc_plastic(n_compositional_fields);
              std::vector<double> visc_effective(n_compositional_fields);

              // Loop over the compositions to compute their viscosities.
              for (unsigned int n=0; n<composition.size(); n++)
                {
                  // Compute the inverse of the diffusion viscosity.
                  visc_diffusion_inverse[n] = 1.0 / diffusion(prefactors_diffusion_fields[n],
                                                              activation_energies_diffusion_fields[n],
                                                              activation_volumes_diffusion_fields[n],
                                                              temperature,
                                                              pressure,
                                                              nu_diffusion_fields[n]);

                  // Compute the inverse of the dislocation viscosity.
                  visc_dislocation_inverse[n] = 1.0 / dislocation(prefactors_dislocation_fields[n],
                                                                  stress_exponents_fields[n],
                                                                  activation_energies_dislocation_fields[n],
                                                                  activation_volumes_dislocation_fields[n],
                                                                  temperature,
                                                                  pressure,
                                                                  strain_rate_dev_inv,
                                                                  nu_dislocation_fields[n]);

                  // Compute the composite viscosity.
                  visc_viscous[n] = 1.0 / (visc_diffusion_inverse[n] + visc_dislocation_inverse[n]);

                  // Compute the plastic viscosity.
                  visc_plastic[n] =  plastic(phis_fields[n],
                                             cohesions_fields[n],
                                             pressure,
                                             strain_rate_dev_inv);

                  // Apply strain rate weakening of plastic viscosity if needed.
                  if (strain_rate_weakening)
                    visc_plastic[n] *= (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1));

                  // Either average the plastic and viscous viscosity harmonically (Eq. (12)),
                  // or take the minimum (Eq. (11)).
                  if (harmonic_plastic_viscous)
                    visc_effective[n] = 1.0 / ((1.0/visc_plastic[n]) + (1.0/visc_viscous[n]));
                  else
                    visc_effective[n] = std::min(visc_viscous[n],visc_plastic[n]);
                }

              // Averaging of the viscosities of each composition
              if (viscosity_averaging == "Harmonic")
                viscosity = harmonic_average(composition, visc_effective);
              else if (viscosity_averaging == "Arithmetic")
                viscosity = arithmetic_average(composition,visc_effective);
              else
                viscosity = geometric_average(composition,visc_effective);

            }

        }
      // No additional compositional fields
      else
        {
          // Compute the inverse of the diffusion viscosity.
          const double visc_diffusion_inverse = 1.0/diffusion(prefactor_diffusion,
                                                              activation_energy_diffusion,
                                                              activation_volume_diffusion,
                                                              temperature,
                                                              pressure,
                                                              nu_diffusion);

          // Compute the inverse of the dislocation viscosity.
          const double visc_dislocation_inverse = 1.0/dislocation(prefactor_dislocation,
                                                                  stress_exponent,
                                                                  activation_energy_dislocation,
                                                                  activation_volume_dislocation,
                                                                  temperature,
                                                                  pressure,
                                                                  strain_rate.norm(),
                                                                  nu_dislocation);

          // Compute the composite viscosity.
          viscosity_viscous = 1.0 / (visc_diffusion_inverse + visc_dislocation_inverse);

          // Compute the plastic viscosity.
          viscosity_plastic = plastic(phi,
                                      C,
                                      pressure,
                                      strain_rate_dev_inv);

          // Apply strain rate weakening of plastic viscosity if needed.
          if (strain_rate_weakening)
            viscosity_plastic *= (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1));

          // Either average the plastic and viscous viscosity harmonically (Eq. (12)),
          // or take the minimum (Eq. (11)).
          if (harmonic_plastic_viscous)
            viscosity = 1.0 / ((1.0/viscosity_plastic) + (1.0/viscosity_viscous));
          else
            viscosity = std::min(viscosity_plastic,viscosity_viscous);
        }

      // If an additional weak zone is specified, apply it.
      if (weak_zone)
        viscosity /= weak_zone_function.value(position);

      // Either take the harmonic average of the viscosity and
      // maximum viscosity (Eq. (14)), or perform a true cut-off (Eq. (13)).
      if (harmonic_max)
        {
          viscosity = 1.0 / ((1.0 / viscosity) + (1.0 / maximum_eta));
          viscosity += minimum_eta;
        }
      else
        viscosity = std::max(std::min(viscosity,maximum_eta),minimum_eta);

      return viscosity;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double> &composition,
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &position) const
    {
      // The effective strain rate
      const SymmetricTensor<2,dim> strain_rate_dev = deviator(strain_rate);
      const double strain_rate_dev_inv = std::sqrt(0.5) * strain_rate_dev.norm();

      // The ratio of plastic to viscous and diffusion to dislocation viscosity.
      double ratio_plastic = 0.0, ratio_composite = 0.0;

      // Multiple compositional fields
      if (composition.size() > 0)
        {
          if (viscosity_averaging == "Max")
            {
              // The index of the compositional field with the highest value.
              const int max_comp = maximum_composition(composition);

              // Compute the inverse of the dislocation viscosity.
              const double visc_dislocation = dislocation(prefactors_dislocation_fields[max_comp],
                                                          stress_exponents_fields[max_comp],
                                                          activation_energies_dislocation_fields[max_comp],
                                                          activation_volumes_dislocation_fields[max_comp],
                                                          temperature,
                                                          pressure,
                                                          strain_rate_dev_inv,
                                                          nu_dislocation_fields[max_comp]);

              // Compute the inverse of the diffusion viscosity.
              const double visc_diffusion = diffusion(prefactors_diffusion_fields[max_comp],
                                                      activation_energies_diffusion_fields[max_comp],
                                                      activation_volumes_diffusion_fields[max_comp],
                                                      temperature,
                                                      pressure,
                                                      nu_diffusion_fields[max_comp]);

              // Compute the plastic viscosity.
              const double visc_plastic = plastic(phis_fields[max_comp],
                                                  cohesions_fields[max_comp],
                                                  pressure,
                                                  strain_rate_dev_inv)
                                          * (strain_rate_weakening ?
                                             (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1)) :
                                             1.0);

              // Compute the ratio of dislocation to diffusion viscosity.
              ratio_composite = visc_dislocation / visc_diffusion;

              // Compute the ratio of composite to plastic viscosity.
              ratio_plastic = (1.0 / ((1.0/visc_dislocation) + (1.0/visc_diffusion))) / visc_plastic;

            }
          // Harmonic, arithmetic and geometric averaging
          else
            {
              std::vector<double> visc_diffusion(n_compositional_fields);
              std::vector<double> visc_dislocation(n_compositional_fields);
              std::vector<double> visc_viscous(n_compositional_fields);
              std::vector<double> visc_plastic(n_compositional_fields);
              for (unsigned int n = 0; n < composition.size(); ++n)
                {
                  // Compute the inverse of the dislocation viscosity.
                  visc_dislocation[n] = dislocation(prefactors_dislocation_fields[n],
                                                    stress_exponents_fields[n],
                                                    activation_energies_dislocation_fields[n],
                                                    activation_volumes_dislocation_fields[n],
                                                    temperature,
                                                    pressure,
                                                    strain_rate_dev_inv,
                                                    nu_dislocation_fields[n]);

                  // Compute the inverse of the dislocation viscosity.
                  visc_diffusion[n] = diffusion(prefactors_diffusion_fields[n],
                                                activation_energies_diffusion_fields[n],
                                                activation_volumes_diffusion_fields[n],
                                                temperature,
                                                pressure,
                                                nu_diffusion_fields[n]);

                  // Compute the composite viscosity.
                  visc_viscous[n] = 1.0 / ((1.0 / visc_dislocation[n]) + (1.0 / visc_diffusion[n]));

                  // Compute the plastic viscosity.
                  visc_plastic[n] = plastic(phis_fields[n],
                                            cohesions_fields[n],
                                            pressure,
                                            strain_rate_dev_inv)
                                    * (strain_rate_weakening ?
                                       (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1)) :
                                       1.0);
                }

              double viscosity_dislocation = 0.0, viscosity_diffusion = 0.0;
              double viscosity_plastic = 0.0, viscosity_viscous = 0.0;

              // Harmonic averaging of the viscosities of each composition
              if (viscosity_averaging == "Harmonic")
                {
                  viscosity_plastic = harmonic_average(composition, visc_plastic);
                  viscosity_dislocation = harmonic_average(composition, visc_dislocation);
                  viscosity_diffusion = harmonic_average(composition, visc_diffusion);
                  viscosity_viscous = harmonic_average(composition, visc_viscous);
                }
              // Arithmetic averaging
              else if (viscosity_averaging == "Arithmetic")
                {
                  viscosity_plastic = arithmetic_average(composition, visc_plastic);
                  viscosity_dislocation = arithmetic_average(composition, visc_dislocation);
                  viscosity_diffusion = arithmetic_average(composition, visc_diffusion);
                  viscosity_viscous = arithmetic_average(composition, visc_viscous);
                }
              // Geometric averaging
              else
                {
                  viscosity_plastic = geometric_average(composition, visc_plastic);
                  viscosity_dislocation = geometric_average(composition, visc_dislocation);
                  viscosity_diffusion = geometric_average(composition, visc_diffusion);
                  viscosity_viscous = geometric_average(composition, visc_viscous);
                }

              // Compute the ratio of composite to plastic viscosity.
              ratio_plastic = viscosity_viscous / viscosity_plastic;

              // Compute the ratio of dislocation to diffusion viscosity.
              ratio_composite = viscosity_dislocation / viscosity_diffusion;
            }
        }
      // No compositional fields
      else
        {
          // Compute the dislocation viscosity.
          const double visc_dislocation = dislocation(prefactor_dislocation,
                                                      stress_exponent,
                                                      activation_energy_dislocation,
                                                      activation_volume_dislocation,
                                                      temperature,
                                                      pressure,
                                                      strain_rate_dev_inv,
                                                      nu_dislocation);

          // Compute the diffusion viscosity.
          const double visc_diffusion = diffusion(prefactor_diffusion,
                                                  activation_energy_diffusion,
                                                  activation_volume_diffusion,
                                                  temperature,
                                                  pressure,
                                                  nu_diffusion);

          // Compute the plastic viscosity.
          const double viscosity_plastic = plastic(phi,
                                                   C,
                                                   pressure,
                                                   strain_rate_dev_inv)
                                           * (strain_rate_weakening ?
                                              (std::max(1.0-(strain_rate_dev_inv/ref_strain_rate),0.1)) :
                                              1.0);

          // Compute the composite viscosity.
          const double viscosity_viscous = 1.0 / ((1.0 / visc_dislocation) + (1.0 / visc_diffusion));

          // Compute the ratio of composite to plastic viscosity.
          ratio_plastic = viscosity_viscous / viscosity_plastic;

          // Compute the ratio of dislocation to diffusion viscosity.
          ratio_composite = visc_dislocation / visc_diffusion;
        }

      double ratio = 0;
      // Dislocation creep is dominant, return 0.
      if (ratio_plastic <= 1.0 && ratio_composite <= 1.0)
        ratio = 0.0;
      // Diffusion creep is dominant, return -1.
      else if (ratio_plastic <= 1.0 && ratio_composite > 1.0)
        ratio = -1.0;
      // Plastic deformation is dominant, return 1.
      else ratio = 1.0;

      return ratio;
    }


    template <int dim>
    double
    Viscoplastic<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &composition,
                   const Point<dim> &) const
    {
      double cp_total=0.0;
      // In case the viscosity uses the infinity norm,
      // also use it here to average the specific heat.
      // Otherwise, use an arithmetic average.
      if (composition.size()>0)
        {
          if (viscosity_averaging == "Max")
            {
              const int max_comp = maximum_composition(composition);
              cp_total = capacities_fields[max_comp];
            }
          else
            cp_total = arithmetic_average(composition, capacities_fields);
        }
      else
        cp_total = reference_specific_heat;

      return cp_total;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &composition,
                          const Point<dim> &) const
    {
      double k_total=0.0;
      // In case the viscosity uses the infinity norm,
      // also use it here to average the thermal conductivity.
      // Otherwise, use an arithmetic average.
      if (composition.size()>0)
        {
          if (viscosity_averaging == "Max")
            {
              const int max_comp = maximum_composition(composition);
              k_total = conductivities_fields[max_comp];
            }
          else
            k_total = arithmetic_average(composition, conductivities_fields);
        }
      else
        k_total = k_value;

      return k_total;
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    Viscoplastic<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &composition,
             const Point<dim> &) const
    {
      double refrho_total=0.0;
      double refT_total=0.0;

      // In case the viscosity uses the infinity norm,
      // also use it here to average the thermal conductivity.
      // Otherwise, use an arithmetic average.
      if (composition.size()>0)
        {
          if (viscosity_averaging=="Max")
            {
              const int max_comp = maximum_composition(composition);
              refrho_total = refdens_fields[max_comp];
              refT_total = reftemps_fields[max_comp];
            }
          else
            {
              refrho_total = arithmetic_average(composition, refdens_fields);
              refT_total   = arithmetic_average(composition, reftemps_fields);
            }
        }
      else
        {
          refrho_total = reference_rho;
          refT_total   = reference_T;
        }

      return (refrho_total * (1.0 - thermal_alpha * (temperature - refT_total)));
    }


    template <int dim>
    double
    Viscoplastic<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &,
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    Viscoplastic<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &,
                     const Point<dim> &) const
    {
      return 0.0;
    }

    template <int dim>
    bool
    Viscoplastic<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (n_compositional_fields != 0)
        {
          return ((dependence & NonlinearDependence::pressure)
                  ||
                  (dependence & NonlinearDependence::strain_rate)
                  ||
                  (dependence & NonlinearDependence::temperature)
                  ||
                  (dependence & NonlinearDependence::compositional_fields));
        }
      else
        {
          return ((dependence & NonlinearDependence::pressure)
                  ||
                  (dependence & NonlinearDependence::strain_rate)
                  ||
                  (dependence & NonlinearDependence::temperature));
        }
    }


    template <int dim>
    bool
    Viscoplastic<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
          &&
          (n_compositional_fields != 0))
        return true;
      else if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
               &&
               (thermal_alpha != 0))
        return true;
      else
        return false;

    }

    template <int dim>
    bool
    Viscoplastic<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Viscoplastic<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
          &&
          (n_compositional_fields != 0))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Viscoplastic<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
          &&
          (n_compositional_fields != 0))
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    Viscoplastic<dim>::
    is_compressible () const
    {
      return false;
    }


    template <int dim>
    void
    Viscoplastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional fields");
      {
        prm.declare_entry ("Number of fields", "0",
                           Patterns::Integer (0),
                           "The number of fields that will be advected along with the flow field, excluding "
                           "velocity, pressure and temperature.");

        // Parameters when compostional fields are used
        prm.declare_entry ("List of phis of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of angles of internal friction equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of conductivities of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of thermal conductivities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of capacities of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of heat capacities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of reftemps of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of reference temperatures equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of refdens of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of reference densities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of initial viscs of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of initial viscosities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of cohesions of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of cohesions equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of prefactors diffusion of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of prefactors of diffusion equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of activation energies diffusion of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of activation energies of diffusion equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of activation volumes diffusion of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of activation volumes of diffusion equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of constant coefficients nu diffusion of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of constant coefficients of diffusion equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of prefactors dislocation of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of prefactors of dislocation equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of constant coefficients nu dislocation of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of constant coefficients of dislocation equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of activation energies dislocation of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of activation energies of dislocation equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of activation volumes dislocation of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of activation volumes of dislocation equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of stress exponents of fields", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of stress exponents equal to the number of "
                           "compositional fields.");
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoplastic model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$ (for normalization). Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Reference viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the reference constant viscosity (for normalization). Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the reference thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the reference specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");

          prm.enter_subsection ("Viscosity");
          {
            prm.declare_entry ("Viscosity Averaging", "Harmonic",
                               Patterns::Anything (),
                               "Averaging of compositional field contributions to viscosity, "
                               "density, specific heat and thermal conductivity.");
            prm.declare_entry ("Harmonic viscous and plastic viscosity averaging", "true",
                               Patterns::Bool (),
                               "Averaging of viscous and plastic viscosity. Can be harmonic or minimum.");
            prm.declare_entry ("Harmonic effective and maximum viscosity averaging", "false",
                               Patterns::Bool (),
                               "Averaging of effective and maximum viscosity. Can be min/max or harmonic "
                               "max viscosity + adding "
                               "the minimum viscosity. The latter option is more smooth but can give "
                               "numerical breakdown.");
            prm.declare_entry ("Initial Viscosity", "5e22",
                               Patterns::Double (0),
                               "The value of the initial viscosity. Units: $kg/m/s$.");
            prm.declare_entry ("Minimum Viscosity", "1e20",
                               Patterns::Double (0),
                               "The value of the minimum viscosity cutoff. Units: $kg/m/s$.");
            prm.declare_entry ("Maximum Viscosity", "1e25",
                               Patterns::Double (0),
                               "The value of the maximum viscosity cutoff. Units: $kg/m/s$.");
            prm.declare_entry ("Use initial strain rate", "false",
                               Patterns::Bool (),
                               "Whether to use an initial strain rate to calculate the viscosity during "
                               "the first nonlinear iteration of the first timestep (true) or to use "
                               "the initial viscosities specified for each compositional field. ");
            prm.declare_entry ("Initial strain rate", "1e-15",
                               Patterns::Double (0),
                               "The strain rate value to use as an initial strain rate to calculate the viscosity during "
                               "the first nonlinear iteration of the first timestep. ");
            prm.declare_entry ("Weak zone", "false",
                               Patterns::Bool (),
                               "Presence of a weak zone.");
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
            prm.declare_entry ("Correct for uniaxial measurements", "false",
                               Patterns::Bool (),
                               "Whether or not to correct dislocation and diffusion prefactors for "
                               "the fact that measurements were performed in a uniaxial setup. "
                               "If true, the prefactors are multiplied with 3^((n+1)/2)*1/2.");
            prm.declare_entry ("Activation energy diffusion", "335e3",
                               Patterns::Double (0),
                               "Activation energy for diffusion creep");
            prm.declare_entry ("Activation volume diffusion", "4.0e-6",
                               Patterns::Double (0),
                               "Activation volume for diffusion creep");
            prm.declare_entry ("Prefactor diffusion", "1.92e-11",
                               Patterns::Double (0),
                               "Prefactor for diffusion creep "
                               "(1e0/prefactor)*exp((activation_energy+activation_volume*pressure)/(R*temperature))");
            prm.declare_entry ("Constant coefficient diffusion", "1.0",
                               Patterns::Double (0),
                               "Constant coefficient for diffusion creep");
            prm.declare_entry ("Activation energy dislocation", "540e3",
                               Patterns::Double (0),
                               "Activation energy for dislocation creep");
            prm.declare_entry ("Activation volume dislocation", "14.0e-6",
                               Patterns::Double (0),
                               "Activation volume for dislocation creep");
            prm.declare_entry ("Prefactor dislocation", "2.42e-10",
                               Patterns::Double (0),
                               "Prefactor for dislocation creep "
                               "(1e0/prefactor)*exp((activation_energy+activation_volume*pressure)/(R*temperature))");
            prm.declare_entry ("Stress exponent", "3.5",
                               Patterns::Double (0),
                               "Stress exponent for dislocation creep");
            prm.declare_entry ("Constant coefficient dislocation", "1.0",
                               Patterns::Double (0),
                               "Constant coefficient for dislocation creep");
            prm.declare_entry ("Angle internal friction", "20",
                               Patterns::Double (0),
                               "Angle of internal friction for plastic creep in absence of compositional fields");
            prm.declare_entry ("Cohesion", "20.0e6",
                               Patterns::Double (0),
                               "Cohesion for plastic creep in absence of compositional fields");
            prm.declare_entry ("Plastic strain rate weakening", "false",
                               Patterns::Bool (),
                               "Whether or not to use strain rate weakening in the calculation of plastic viscosity");
            prm.declare_entry ("Reference strain rate", "1e-15",
                               Patterns::Double (0),
                               "The reference strain rate used in strain rate weakening");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Viscoplastic<dim>::parse_parameters (ParameterHandler &prm)
    {


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoplastic model");
        {

          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          reference_eta              = prm.get_double ("Reference viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");


          prm.enter_subsection ("Viscosity");
          {
            viscosity_averaging      = prm.get ("Viscosity Averaging");
            AssertThrow(viscosity_averaging == "Max" || viscosity_averaging == "Harmonic" ||
                        viscosity_averaging == "Geometric" || viscosity_averaging == "Arithmetic",
                        ExcMessage("Invalid input parameter file: This type of viscosity averaging is not implemented"));
            harmonic_plastic_viscous = prm.get_bool ("Harmonic viscous and plastic viscosity averaging");
            harmonic_max      = prm.get_bool ("Harmonic effective and maximum viscosity averaging");
            initial_eta       = prm.get_double ("Initial Viscosity");
            minimum_eta       = prm.get_double ("Minimum Viscosity");
            maximum_eta       = prm.get_double ("Maximum Viscosity");
            use_initial_strain_rate = prm.get_bool ("Use initial strain rate");
            initial_strain_rate = prm.get_double ("Initial strain rate");
            weak_zone         = prm.get_bool   ("Weak zone");
            weak_zone_function.parse_parameters (prm);
            correct_uniaxiality = prm.get_bool ("Correct for uniaxial measurements");

            activation_energy_diffusion   = prm.get_double ("Activation energy diffusion");
            activation_volume_diffusion   = prm.get_double ("Activation volume diffusion");
            prefactor_diffusion           = prm.get_double ("Prefactor diffusion") * (correct_uniaxiality ? (3. / 2.) : 1.);
            nu_diffusion                  = prm.get_double ("Constant coefficient diffusion");
            activation_energy_dislocation = prm.get_double ("Activation energy dislocation");
            activation_volume_dislocation = prm.get_double ("Activation volume dislocation");
            prefactor_dislocation         = prm.get_double ("Prefactor dislocation") * (correct_uniaxiality ? (std::pow(3,(stress_exponent+1.0)/2.0)*0.5) : 1.);
            nu_dislocation                = prm.get_double ("Constant coefficient dislocation");
            stress_exponent               = prm.get_double ("Stress exponent");
            phi                           = prm.get_double ("Angle internal friction");
            C                             = prm.get_double ("Cohesion");
            strain_rate_weakening         = prm.get_bool ("Plastic strain rate weakening");
            ref_strain_rate               = prm.get_double ("Reference strain rate");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.leave_subsection();

      prm.enter_subsection ("Compositional fields");
      {
        n_compositional_fields = prm.get_integer ("Number of fields");

        if (n_compositional_fields > 0)
          {

            // Parameters needed for all rheologies
            const std::vector<double> n_conductivities_fields = Utilities::string_to_double
                                                                (Utilities::split_string_list(prm.get ("List of conductivities of fields")));
            conductivities_fields = std::vector<double> (n_conductivities_fields.begin(),
                                                         n_conductivities_fields.end());
            AssertThrow (conductivities_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of conductivities of fields"));


            const std::vector<double> n_capacities_fields = Utilities::string_to_double
                                                            (Utilities::split_string_list(prm.get ("List of capacities of fields")));
            capacities_fields = std::vector<double> (n_capacities_fields.begin(),
                                                     n_capacities_fields.end());
            AssertThrow (capacities_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of capacities of fields"));

            const std::vector<double> n_reftemps_fields = Utilities::string_to_double
                                                          (Utilities::split_string_list(prm.get ("List of reftemps of fields")));
            reftemps_fields = std::vector<double> (n_reftemps_fields.begin(),
                                                   n_reftemps_fields.end());
            AssertThrow (reftemps_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of reftemps of fields"));


            const std::vector<double> n_refdens_fields = Utilities::string_to_double
                                                         (Utilities::split_string_list(prm.get ("List of refdens of fields")));
            refdens_fields = std::vector<double> (n_refdens_fields.begin(),
                                                  n_refdens_fields.end());
            AssertThrow (refdens_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of refdens of fields"));

            const std::vector<double> n_init_eta_fields = Utilities::string_to_double
                                                          (Utilities::split_string_list(prm.get ("List of initial viscs of fields")));
            init_eta_fields = std::vector<double> (n_init_eta_fields.begin(),
                                                   n_init_eta_fields.end());
            AssertThrow (init_eta_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of initial viscs of fields"));

            // Plastic parameters
            const std::vector<double> n_cohesions_fields = Utilities::string_to_double
                                                           (Utilities::split_string_list(prm.get ("List of cohesions of fields")));
            cohesions_fields = std::vector<double> (n_cohesions_fields.begin(),
                                                    n_cohesions_fields.end());
            AssertThrow (cohesions_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of cohesions of fields"));

            const std::vector<double> n_phis_fields = Utilities::string_to_double
                                                      (Utilities::split_string_list(prm.get ("List of phis of fields")));
            // Convert friction angles to degrees
            for (unsigned int i=0; i<n_compositional_fields; i++)
              phis_fields.push_back(n_phis_fields[i]*numbers::PI/180.0);

            AssertThrow (phis_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of phis of fields"));

            // Dislocation parameters
            const std::vector<double> n_stressexponent_fields = Utilities::string_to_double
                                                                (Utilities::split_string_list(prm.get ("List of stress exponents of fields")));
            stress_exponents_fields = std::vector<double> (n_stressexponent_fields.begin(),
                                                           n_stressexponent_fields.end());
            AssertThrow (stress_exponents_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of stress exponents of fields"));

            const std::vector<double> n_prefactordisl_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of prefactors dislocation of fields")));
            // Correction of measurements for uniaxiality
            for (unsigned int i=0; i<n_compositional_fields; i++)
              {
                if (correct_uniaxiality)
                  prefactors_dislocation_fields.push_back(n_prefactordisl_fields[i]*std::pow(3,(n_stressexponent_fields[i]+1.0)/2.0)*0.5);
                else
                  prefactors_dislocation_fields.push_back(n_prefactordisl_fields[i]);
              }
            AssertThrow (prefactors_dislocation_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of prefactors dislocation of fields"));

            const std::vector<double> n_actenergydisl_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of activation energies dislocation of fields")));
            activation_energies_dislocation_fields = std::vector<double> (n_actenergydisl_fields.begin(),
                                                                          n_actenergydisl_fields.end());
            AssertThrow (activation_energies_dislocation_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of activation energies dislocation of fields"));

            const std::vector<double> n_actvolumedisl_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of activation volumes dislocation of fields")));
            activation_volumes_dislocation_fields = std::vector<double> (n_actvolumedisl_fields.begin(),
                                                                         n_actvolumedisl_fields.end());
            AssertThrow (activation_volumes_dislocation_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of activation volumes dislocation of fields"));

            const std::vector<double> n_nu_fields = Utilities::string_to_double
                                                    (Utilities::split_string_list(prm.get ("List of constant coefficients nu dislocation of fields")));
            nu_dislocation_fields = std::vector<double> (n_nu_fields.begin(),
                                                         n_nu_fields.end());
            AssertThrow (nu_dislocation_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of constant coefficients nu dislocation of fields"));

            // Diffusion parameters
            const std::vector<double> n_prefactordiff_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of prefactors diffusion of fields")));
            // Correction for uniaxial measurements
            for (unsigned int i=0; i<n_compositional_fields; i++)
              {
                if (correct_uniaxiality)
                  prefactors_diffusion_fields.push_back(n_prefactordiff_fields[i]*3.0*0.5);
                else
                  prefactors_diffusion_fields.push_back(n_prefactordiff_fields[i]);
              }
            AssertThrow (prefactors_diffusion_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of prefactors diffusion of fields"));

            const std::vector<double> n_actenergydiff_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of activation energies diffusion of fields")));
            activation_energies_diffusion_fields = std::vector<double> (n_actenergydiff_fields.begin(),
                                                                        n_actenergydiff_fields.end());
            AssertThrow (activation_energies_diffusion_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of activation energies diffusion of fields"));

            const std::vector<double> n_actvolumediff_fields = Utilities::string_to_double
                                                               (Utilities::split_string_list(prm.get ("List of activation volumes diffusion of fields")));
            activation_volumes_diffusion_fields = std::vector<double> (n_actvolumediff_fields.begin(),
                                                                       n_actvolumediff_fields.end());
            AssertThrow (activation_volumes_diffusion_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of activation volumes diffusion of fields"));

            const std::vector<double> n_nu_diff_fields = Utilities::string_to_double
                                                         (Utilities::split_string_list(prm.get ("List of constant coefficients nu diffusion of fields")));
            nu_diffusion_fields = std::vector<double> (n_nu_diff_fields.begin(),
                                                       n_nu_diff_fields.end());
            AssertThrow (nu_diffusion_fields.size() == n_compositional_fields,
                         ExcMessage("Invalid input parameter file: Wrong number of entries in List of constant coefficients nu diffusion of fields"));

          }
      }
      prm.leave_subsection();

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Viscoplastic,
                                   "viscoplastic",
                                   "A material model that allows for nonlinear rheologies for none, 1 or more compositions."
                                   "Viscosity can depend on temperature, pressure, strain rate and composition and "
                                   "constitutes a combination of Drucker-Prager plasticity, dislocation creep and diffusion creep. "
                                   "Averaging of the viscosity contributions of the compositional fields occurs through "
                                   "harmonic, geometric, arithmetic or infinity norm averaging."
                                   "Density depends on temperature, reference temperature and reference density, which depend on composition. "
                                   "Specific heat and thermal conductivity can depend on composition, but are otherwise constant. "
                                   "Note that material parameters per composition field are specified in a different subsection than "
                                   "the material parameters used when no fields are present.")
  }
}
