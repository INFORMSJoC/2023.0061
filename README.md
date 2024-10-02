[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An iterative exact algorithm over a time-expanded network for the transportation of biomedical samples

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper An iterative exact algorithm over a time-expanded network for the transportation of biomedical samples by Daniel M. Ocampo-Giraldo, Ana M. Anaya-Arenas and Claudio Contardo. The upstream repository (including a potentially more up-to-date version of this code) can be found at https://github.com/dmocampog/BSTPDDD.jl. This snapshot has been built from [this commit](https://github.com/dmocampog/BSTPDDD.jl/commit/c42abe2ae1b0ba955300cbba5397d8e5fc4eb819)


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0061

https://doi.org/10.1287/ijoc.2023.0061.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{bstp2024ijoc,
  author =        {D. M. Ocampo-Giraldo and A. M. Anaya-Arenas and Claudio Contardo},
  publisher =     {INFORMS Journal on Computing},
  title =         {An iterative exact algorithm over a time-expanded network for the transportation of biomedical samples},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0061.cd},
  url =           {https://github.com/INFORMSJoC/2023.0061},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0061},
}
```

## Description

This repository provides the source code, data files and detailed results as reported in the article.

The main folders are 'data', 'results', 'src' and 'test'.
- '[inst](inst)': This folder includes instances of the BSTP.
- '[src](src)': The source code.
- '[test](test)': Tests.


## Installation

To install this Julia package, simply execute from a Julia Pkg REPL (by pressing `]` within a regular Julia REPL) the following:
```julia
(@v1.10) pkg> add https://github.com/INFORMSJoC/2023.0061
```

## Testing

You can test this method by executing from a Julia Pkg REPL (by pressing `]` within a regular Julia REPL) the following:
```julia
(@v1.10) pkg> test BSTPDDD
```

## Basic Usage

To solve instance `I01` (which must be present inside the folder `inst`) for $\Delta_k = 15$ and $\Delta_t = 5$ please execute
```julia
julia> using BSTPDDD
julia> solution = BSTPDDD.ddd_bstp("I01",15,5)

# returns:
#    Number of SCCs
#    Number of Transportation Request
#    Solution routes
#    Uperbund value
#    Lowerbound value
#    Optimality gap
#    CPU time
#    Number of iterations
#    Termination criteria

# solution routes structure
#  Commodities group id picked-up by the route
#  SCCs visited by the route
#  Arriving time to first SCCs in the route
#  Departing time form last SCCs in the route
#  Total route's slack time
```
