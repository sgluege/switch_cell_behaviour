# switch_cell_behaviour

Switch behaviour of a cell based on its type. Here a single biomodule (`GeneralModule`) is used, that runs the simulation 
step based on the `cell_type`, either `type1behaviour()` or `type2behaviour()` is called in the modules' `Run()`routine.

Simulation starts with `cell_type=1` and is switched to `cell_type=2` when the  intracellular substance quantity `counter_sub_quant` falls below a predefind threshold.

-----------------
##### BioDynaMo
```
$ biodynamo --version
v0.1.0-371-g7e4f5bd5

```