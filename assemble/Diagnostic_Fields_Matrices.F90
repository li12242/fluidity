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

module diagnostic_fields_matrices
  !!< A module to link to diagnostic variable calculations.

  use global_parameters, only:FIELD_NAME_LEN 
  use fields
  use sparse_matrices_fields
  use field_derivatives
  use state_module
  use futils
  use fetools
  use spud
  use parallel_fields, only: zero_non_owned
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use divergence_matrix_cg, only: assemble_divergence_matrix_cg, assemble_compressible_divergence_matrix_cg
  use gradient_matrix_cg, only: assemble_gradient_matrix_cg
  use parallel_tools
  use sparsity_patterns, only: make_sparsity
  use state_fields_module
  use solvers
  use transform_elements
  use sparse_tools
  use sparsity_patterns_meshes
  use state_matrices_module
  use halos

  implicit none

  private
  
  public :: calculate_divergence_cv, calculate_divergence_fe, &
            calculate_div_t_cv, calculate_div_t_fe, &
            calculate_grad_fe, calculate_sum_velocity_divergence, &
            calculate_compressible_continuity_residual

contains

  subroutine calculate_divergence_cv(state, div)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: div

      type(vector_field), pointer :: field, x

      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(scalar_field), pointer :: cvmass
      type(scalar_field) :: inverse_cvmass
      type(scalar_field) :: ctfield, ct_rhs

      call get_option(trim(div%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_vector_field(state, trim(field_name))
      x=>extract_vector_field(state, "Coordinate")

      call allocate(ctfield, div%mesh, name="CTField")
      call allocate(ct_rhs, div%mesh, name="CTRHS")

      cvmass => get_cv_mass(state, div%mesh)

      CT_m => get_divergence_matrix_cv(state, test_mesh=div%mesh, field=field, div_rhs=ct_rhs)

      call mult(ctfield, CT_m, field)
      call addto(ctfield, ct_rhs, -1.0)

      call allocate(inverse_cvmass, cvmass%mesh, "InverseCVMass")
      call invert(cvmass, inverse_cvmass)
      call set(div, ctfield)
      call scale(div, inverse_cvmass)
      
      call deallocate(inverse_cvmass)
      call deallocate(ctfield)
      call deallocate(ct_rhs)

      call halo_update(div)

  end subroutine calculate_divergence_cv

  subroutine calculate_div_t_cv(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(scalar_field), pointer :: lumped_mass
      type(scalar_field) :: inverse_lumped_mass
      type(csr_matrix), pointer :: mass, inverse_mass
      type(vector_field) :: cfield
      type(scalar_field) :: mag, inverse_mag

      logical :: lump_mass, dg, normalise

      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)
      
      dg = (continuity(grad)<0)
      lump_mass = have_option(trim(grad%option_path)//&
                  &"/diagnostic/lump_mass_matrix")
      normalise = have_option(trim(grad%option_path)//&
                  &"/diagnostic/normalise")
      
      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      CT_m => get_divergence_matrix_cv(state, test_mesh=field%mesh, field=grad, exclude_boundaries=.true.)

      call mult_T(cfield, CT_m, field)
      call scale(cfield, -1.0)
      
      call zero(grad)
      
      if(lump_mass) then
        
        lumped_mass => get_lumped_mass(state, grad%mesh)
        call allocate(inverse_lumped_mass, lumped_mass%mesh, "InverseLumpedMass")
        
        call invert(lumped_mass, inverse_lumped_mass)
        call set(grad, cfield)
        call scale(grad, inverse_lumped_mass)
       
        call deallocate(inverse_lumped_mass)
      else if(dg) then
        inverse_mass => get_dg_inverse_mass(state, grad%mesh)
        call mult(grad, inverse_mass, cfield)
      else
        mass => get_mass_matrix(state, grad%mesh)
        call petsc_solve(grad, mass, cfield)
      end if

      if(normalise) then
        mag = magnitude(grad)
        call allocate(inverse_mag, mag%mesh, "InverseMagnitude")
        
        call invert(mag, inverse_mag, tolerance=epsilon(0.0))
        call scale(grad, inverse_mag)
        
        call deallocate(inverse_mag)
        call deallocate(mag)
      end if
      
      call deallocate(cfield)

  end subroutine calculate_div_t_cv

  subroutine calculate_divergence_fe(state, div)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: div

      type(vector_field), pointer :: field, x

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(scalar_field) :: ctfield, ct_rhs

      call get_option(trim(div%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_vector_field(state, trim(field_name))
      x=>extract_vector_field(state, "Coordinate")

      call allocate(ctfield, div%mesh, name="CTField")

      divergence_sparsity=make_sparsity(div%mesh, field%mesh, "DivergenceSparsity")
      allocate(CT_m)
      call allocate(CT_m, divergence_sparsity, (/1, field%dim/), name="DivergenceMatrix" )
      call allocate(ct_rhs, div%mesh, name="CTRHS")


      mass_sparsity=make_sparsity(div%mesh, div%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")
      call zero(mass)

      call assemble_divergence_matrix_cg(CT_m, state, ct_rhs=ct_rhs, &
                                         test_mesh=div%mesh, field=field, &
                                         option_path=div%option_path, div_mass=mass)

      call mult(ctfield, CT_m, field)
      call addto(ctfield, ct_rhs, -1.0)

      call zero(div)
      call petsc_solve(div, mass, ctfield)

      call deallocate(CT_m)
      deallocate(CT_m)
      call deallocate(ct_rhs)
      call deallocate(ctfield)
      call deallocate(divergence_sparsity)
      call deallocate(mass_sparsity)
      call deallocate(mass)

  end subroutine calculate_divergence_fe

  subroutine calculate_div_t_fe(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(vector_field) :: cfield


      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      ! Sparsity of C^T - the transpose of the gradient operator.
      divergence_sparsity=make_sparsity(field%mesh, grad%mesh, "DivergenceSparsity")
      allocate(CT_m)
      call allocate(CT_m, divergence_sparsity, (/1, grad%dim/), name="DivergenceMatrix" )

      mass_sparsity=make_sparsity(grad%mesh, grad%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")

      call assemble_divergence_matrix_cg(CT_m, state, &
                                         test_mesh=field%mesh, field=grad, &
                                         grad_mass=mass)

      call mult_T(cfield, CT_m, field)
      call scale(cfield, -1.0)

      call zero(grad)
      call petsc_solve(grad, mass, cfield)

      call deallocate(divergence_sparsity)
      call deallocate(CT_m)
      deallocate(CT_m)
      call deallocate(mass_sparsity)
      call deallocate(mass)
      call deallocate(cfield)

  end subroutine calculate_div_t_fe

  subroutine calculate_grad_fe(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(csr_sparsity) :: gradient_sparsity
      type(block_csr_matrix) :: C_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(vector_field) :: cfield


      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      ! Sparsity of C^T - the transpose of the gradient operator.
      gradient_sparsity=make_sparsity(grad%mesh, field%mesh, "GradientSparsity")
      call allocate(C_m, gradient_sparsity, (/grad%dim, 1/), name="GradientMatrix" )

      mass_sparsity=make_sparsity(grad%mesh, grad%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")

      call assemble_gradient_matrix_cg(C_m, state, &
                                         test_mesh=grad%mesh, field=field, &
                                         option_path=trim(grad%option_path), &
                                         grad_mass=mass)

      call mult(cfield, C_m, field)

      call zero(grad)
      call petsc_solve(grad, mass, cfield)

      call deallocate(gradient_sparsity)
      call deallocate(C_m)
      call deallocate(mass_sparsity)
      call deallocate(mass)
      call deallocate(cfield)

  end subroutine calculate_grad_fe


  subroutine calculate_sum_velocity_divergence(state, sum_velocity_divergence)
      !!< Calculates \sum{div(vfrac*u)}, where we sum over each prognostic velocity 
      !!< field (i.e. each phase). Used in multiphase flow simulations.

      type(state_type), dimension(:), intent(inout) :: state
      type(scalar_field), pointer :: sum_velocity_divergence

      ! Local variables
      type(vector_field), pointer :: u, x
      integer :: i, stat

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix), pointer :: ct_m

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(scalar_field) :: ctfield, ct_rhs, temp
      type(scalar_field), pointer :: cv_mass
      
      logical :: test_with_cv_dual

      ewrite(1,*) 'Entering calculate_sum_velocity_divergence'
         
      ! Are we testing the divergence with the CV dual mesh
      test_with_cv_dual = have_option(trim(sum_velocity_divergence%option_path)//'/diagnostic/test_with_cv_dual')   

      ! Allocate memory for matrices and sparsity patterns
      call allocate(ctfield, sum_velocity_divergence%mesh, name="CTField")
      call zero(ctfield)
      call allocate(temp, sum_velocity_divergence%mesh, name="Temp")

      ! Require a mass matrix type if not CV tested equation.
      if (.not. test_with_cv_dual) then
         mass_sparsity=make_sparsity(sum_velocity_divergence%mesh, sum_velocity_divergence%mesh, "MassSparsity")
         call allocate(mass, mass_sparsity, name="MassMatrix")
         call zero(mass)      
      end if 
      
      ! Sum up over the div's
      do i = 1, size(state)
         u => extract_vector_field(state(i), "Velocity", stat)

         ! If there's no velocity then cycle
         if(stat/=0) cycle
         ! If this is an aliased velocity then cycle
         if(aliased(u)) cycle
         ! If the velocity isn't prognostic then cycle
         if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

         ! If velocity field is prognostic, begin calculations below
         x => extract_vector_field(state(i), "Coordinate")

         ! Allocate sparsity patterns, C^T matrix and C^T RHS for current state
         divergence_sparsity=make_sparsity(sum_velocity_divergence%mesh, u%mesh, "DivergenceSparsity")

         allocate(ct_m)
         call allocate(ct_m, divergence_sparsity, (/1, u%dim/), name="DivergenceMatrix" )
         call allocate(ct_rhs, sum_velocity_divergence%mesh, name="CTRHS")

         ! Reassemble C^T matrix here
         if (test_with_cv_dual) then
            call assemble_divergence_matrix_cv(ct_m, state(i), ct_rhs=ct_rhs, &
                              test_mesh=sum_velocity_divergence%mesh, field=u)
         else
            if(i==1) then ! Construct the mass matrix (just do this once)
               call assemble_divergence_matrix_cg(ct_m, state(i), ct_rhs=ct_rhs, &
                              test_mesh=sum_velocity_divergence%mesh, field=u, &
                              option_path=sum_velocity_divergence%option_path, div_mass=mass)
            else
               call assemble_divergence_matrix_cg(ct_m, state(i), ct_rhs=ct_rhs, &
                              test_mesh=sum_velocity_divergence%mesh, field=u, &
                              option_path=sum_velocity_divergence%option_path)
            end if   

         end if

         ! Construct the linear system of equations
         call zero(temp)
         call mult(temp, ct_m, u)
         call addto(temp, ct_rhs, -1.0)

         ! Now add it to the sum
         call addto(ctfield, temp)

         call deallocate(ct_m)
         deallocate(ct_m)
         call deallocate(ct_rhs)
         call deallocate(divergence_sparsity)

      end do

      ! Solve for sum_velocity_divergence
      ! ( = \sum{div(vfrac*u)} for incompressible multiphase flows )            
      if (test_with_cv_dual) then
         ! get the cv mass matrix
         cv_mass => get_cv_mass(state(1), sum_velocity_divergence%mesh)
         call set(sum_velocity_divergence, cv_mass)
         call invert(sum_velocity_divergence)
         call scale(sum_velocity_divergence, ctfield)
      else
         call zero(sum_velocity_divergence)
         call petsc_solve(sum_velocity_divergence, mass, ctfield)
      end if 
      
      ! Deallocate memory
      call deallocate(ctfield)
      call deallocate(temp)
      
      if (.not. test_with_cv_dual) then
         call deallocate(mass_sparsity)
         call deallocate(mass)      
      end if 

      ewrite(1,*) 'Exiting calculate_sum_velocity_divergence'
         
  end subroutine calculate_sum_velocity_divergence

  
  subroutine calculate_compressible_continuity_residual(state, compressible_continuity_residual)
      !!< Calculates the residual of the continity equation used in compressible multiphase flow simulations:
      !!< vfrac_c*d(rho_c)/dt + div(rho_c*vfrac_c*u_c) + \sum_i{ rho_c*div(vfrac_i*u_i) }
      
      type(state_type), dimension(:), intent(inout) :: state
      type(scalar_field), pointer :: compressible_continuity_residual, density, olddensity, vfrac
      type(scalar_field) :: drhodt

      ! Local variables
      type(vector_field), pointer :: u, x
      integer :: i, stat, ele
      logical :: dg

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix), pointer :: ct_m

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(scalar_field) :: ctfield, ct_rhs, temp
      
      type(element_type) :: test_function
      type(element_type), pointer :: compressible_continuity_residual_shape
      integer, dimension(:), pointer :: compressible_continuity_residual_nodes
      real, dimension(:), allocatable :: detwei
      real, dimension(:), allocatable :: drhodt_addto
      
      real :: dt

      ewrite(1,*) 'Entering calculate_compressible_continuity_residual'

      ! Allocate memory for matrices and sparsity patterns
      call allocate(ctfield, compressible_continuity_residual%mesh, name="CTField")
      call zero(ctfield)
      call allocate(temp, compressible_continuity_residual%mesh, name="Temp")
            
      mass_sparsity=make_sparsity(compressible_continuity_residual%mesh, compressible_continuity_residual%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")
      call zero(mass)
      
      call get_option("/timestepping/timestep", dt)
      
      ! Sum up over the div's
      do i = 1, size(state)
         u => extract_vector_field(state(i), "Velocity", stat)

         ! If there's no velocity then cycle
         if(stat/=0) cycle
         ! If this is an aliased velocity then cycle
         if(aliased(u)) cycle
         ! If the velocity isn't prognostic then cycle
         if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

         ! If velocity field is prognostic, begin calculations below
         x => extract_vector_field(state(i), "Coordinate")

         ! Allocate sparsity patterns, C^T matrix and C^T RHS for current state
         divergence_sparsity=make_sparsity(compressible_continuity_residual%mesh, u%mesh, "DivergenceSparsity")
         allocate(ct_m)
         call allocate(ct_m, divergence_sparsity, (/1, u%dim/), name="DivergenceMatrix")
         call allocate(ct_rhs, compressible_continuity_residual%mesh, name="CTRHS")

         ! Reassemble C^T matrix here
         if(i==1) then ! Construct the mass matrix (just do this once)
            call assemble_compressible_divergence_matrix_cg(ct_m, state, i, ct_rhs, div_mass=mass)
         else
            call assemble_compressible_divergence_matrix_cg(ct_m, state, i, ct_rhs)
         end if   

         if(have_option("/material_phase::"//trim(state(i)%name)//"/equation_of_state/compressible")) then
            ! Get the time derivative term for the compressible phase's density, vfrac_c * d(rho_c)/dt

            call allocate(drhodt, compressible_continuity_residual%mesh, name="drhodt")
                  
            density => extract_scalar_field(state(i), "Density", stat)
            olddensity => extract_scalar_field(state(i), "OldDensity", stat)
            vfrac => extract_scalar_field(state(i), "PhaseVolumeFraction")
            
            ! Assumes Density and OldDensity are on the same mesh as Pressure,
            ! as it should be according to the manual.
            call zero(drhodt)
            
            dg = continuity(compressible_continuity_residual) < 0
            
            element_loop: do ele = 1, element_count(compressible_continuity_residual)

               if(.not.dg .or. (dg .and. element_owned(compressible_continuity_residual,ele))) then
               
                  allocate(detwei(ele_ngi(compressible_continuity_residual, ele)))
                  allocate(drhodt_addto(ele_loc(compressible_continuity_residual, ele)))
               
                  compressible_continuity_residual_nodes => ele_nodes(compressible_continuity_residual, ele)
                  compressible_continuity_residual_shape => ele_shape(compressible_continuity_residual, ele)
                  test_function = compressible_continuity_residual_shape
                  
                  call transform_to_physical(x, ele, detwei=detwei)
                  
                  drhodt_addto = shape_rhs(test_function, detwei*(ele_val_at_quad(vfrac,ele)*(ele_val_at_quad(density,ele) - ele_val_at_quad(olddensity,ele))/dt))
                  
                  call addto(drhodt, compressible_continuity_residual_nodes, drhodt_addto)
                  
                  deallocate(detwei)
                  deallocate(drhodt_addto)
               end if

            end do element_loop
            
            call addto(ctfield, drhodt)    
            
            call deallocate(drhodt)
            
         end if

         ! Construct the linear system of equations
         call zero(temp)
         call mult(temp, ct_m, u)
         call addto(temp, ct_rhs, -1.0)

         ! Now add it to the sum
         call addto(ctfield, temp)

         call deallocate(ct_m)
         deallocate(ct_m)
         call deallocate(ct_rhs)
         call deallocate(divergence_sparsity)

      end do

      ! Solve for compressible_continuity_residual      
      call zero(compressible_continuity_residual)
      call petsc_solve(compressible_continuity_residual, mass, ctfield)
      ewrite_minmax(compressible_continuity_residual)
         
      ! Deallocate memory
      call deallocate(ctfield)
      call deallocate(temp)
      call deallocate(mass_sparsity)
      call deallocate(mass)

      ewrite(1,*) 'Exiting calculate_compressible_continuity_residual'
         
  end subroutine calculate_compressible_continuity_residual
  
  
end module diagnostic_fields_matrices
