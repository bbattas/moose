##############################################################################
# File: simple_input.i
# File Location: /
# Created Date: Sunday April 21st 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Sunday April 21st 2024
# Modified By: Brandon Battas
# -----
# Description:
#
#
#
#
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 1
[]
[Materials]
  [testmat]
    type =
  []
[]
[Executioner]
  type = Transient
  num_steps = 1
[]
[Problem]
  solve = false
[]
