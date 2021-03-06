!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module physics_from_options
use fields
use state_module
use spud
implicit none

private

public get_fs_reference_density_from_options, convert_free_surface_to_pressure

contains

  subroutine get_fs_reference_density_from_options(rho0, phase_path)
    ! The density returned by this routine defines the relation between
    ! free surface elevation and barotropic pressure: p = \rho0 g \eta
    ! For equation types Boussinesq and ShallowWater this is simply 1.0
    ! (corresponding to the constant in front of the du/dt term in the mom. eqn.)
    ! For LinearMomentum, we return the reference density. By default it is indeed
    ! assumed that the density near the f.s. is equal to the reference density and 
    ! thus p = \rho0 g \eta. Density variations are taken into account in assemble/Free_Surface.F90
    ! when using the variable_density option under the fs. bc.

    real, intent(out) :: rho0
    character(len=*), intent(in), optional :: phase_path

    if (have_option(trim(phase_path)//'/vector_field::Velocity/prognostic/equation::ShallowWater') .or. &
        have_option(trim(phase_path)//'/vector_field::Velocity/prognostic/equation::Boussinesq')) then

      ! with Boussinesq \rho0 is already divided out of the mom. eqn. (and therefore out of Pressure)
      ! same goes for ShallowWater
      rho0 = 1.0

    else if (have_option(trim(phase_path)//'/vector_field::Velocity/prognostic/equation::LinearMomentum')) then
      if (have_option(trim(phase_path)//'/equation_of_state/fluids/linear/')) then

        ! for other equation types the pressure is really the pressure so we do need to
        ! divide by the real reference density
        call get_option(trim(phase_path)//'/equation_of_state/fluids/linear/reference_density', rho0)
      else if (have_option(trim(phase_path)//'/equation_of_state/fluids/ocean_pade_approximation/')) then
        ! The reference density is hard-coded to be 1.0 for the Ocean Pade Approximation
        rho0=1.0
      else
        ewrite(-1,*) "Unless using Boussinesq Velocity, you must specify a"
        ewrite(-1,*) "linear or pade equation of state for the free surface."
        FLExit("Error retrieving reference density from options.")
      end if
    else if (.not. have_option(trim(phase_path)//'vector_field::Velocity/prognostic')) then
      ! this should have been options checked somewhere else, so flabort
      FLAbort("Need a prognostic velocity to retrieve free surface reference density.")
    else
      ! this must mean it doesn't have Velocity equation option, or an unknown equation type
      ! both are not allowed in the schema, so flabort
      FLAbort("Unknown equation type for Velocity")
    endif

  end subroutine get_fs_reference_density_from_options

  subroutine convert_free_surface_to_pressure(field, phase_path)
    type(scalar_field), intent(inout) :: field
    character(len=*), intent(in) :: phase_path

    real :: rho0, gravity_magnitude

    call get_fs_reference_density_from_options(rho0, phase_path)

    if(have_option("/physical_parameters/gravity")) then
       call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
       ! Note this is rescale is limited to the region this condition is applied
       call scale(field, gravity_magnitude * rho0)
    else
       !FLExit("Specifying a free surface initial condition requires gravity to be defined.")
       ewrite(-1,*) "Converting free-surface height to a pressure field requires gravity to be defined."
       FLExit("Have you specified a free-surface initial condition without defining gravity?")
    end if

  end subroutine convert_free_surface_to_pressure

end module physics_from_options
