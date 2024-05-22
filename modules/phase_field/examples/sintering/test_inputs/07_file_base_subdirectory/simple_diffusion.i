##############################################################################
# File: simple_diffusion.i
# File Location: /examples/sintering/test_inputs/05_global_variable_issue
# Created Date: Tuesday March 19th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday May 22nd 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Trying to replicate the issue where all the tests i ran using file_base
#   with a name/name value created those folders also in the base moose folder
#   with nothing in them?  Cant seem to replicate it now either way. Maybe
#   it fixed when i just recompiled earlier?
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'hypre'
[]

# [Outputs]
#   # exodus = true
#   nemesis = true
#   file_base = subdirTest/subdirTest
# []

[Outputs]
  perf_graph = false
  csv = true
  exodus = false
  checkpoint = false
  file_base = subdirTest3/subdirTest
  [nemesis]
    type = Nemesis
  []
  print_linear_residuals = false
[]
