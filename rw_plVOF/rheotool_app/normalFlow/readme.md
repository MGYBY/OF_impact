Normal flow simulation for highly non-Newtonian fluids (*n* = 0.20). The fully-developed normal flow is formed at the end of the channel.

2D sim configuration. It may be more reasonable to modify rheoInterFoam to be compatible with the periodic BC (follow [interCyclicFoam](https://github.com/MGYBY/OF_impact/blob/main/interCyclicFoam.zip)).

Normal flow result here agrees with the analytical solution. The possible reasons why built-in rheology model of OF gave unreasonable results:
* Formulation based on kinematic viscosity instead of dynamic viscosity?
* Unbounded alpha
* Rheology model only designed for single-phase flow?
