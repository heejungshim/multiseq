# multiseq: multiscale Poisson process approaches for differential or association analysis of high-throughput sequencing data

This R package implements multiscale Poisson process approaches for differential or association analysis of high-throughput sequencing data. Key features that distinguish the multiseq from typical differential or association analysis are to 1) better exploit high-resolution information in the high-throughput sequencing data, and 2) directly model the count nature of the data. See [Shim et al. (2021)][multiseq-arxiv] for the details of the motivations, approaches, and comparison with other methods. 

If you find a bug, please post an [issue][issues].

## License

Copyright (c) 2021-2023, Heejung Shim, Zhengrong Xing, Ester Pantaleo, and Matthew Stephens.

All source code and software in this repository is free software; you
can redistribute it and/or modify it under the terms of the
[GNU General Public License][gpl] as published by the
[Free Software Foundation][fsf]; either version 3 of the License, or
(at your option) any later version. See the [LICENSE](LICENSE) file
for the full text of the license.

## Citing this work

If you find that this R package useful for your work, please cite our
paper:

> Heejung Shim, Zhengrong Xing, Ester Pantaleo, Francesca Luca, Roger 
> Pique-Regi, and Matthew Stephens (2021). *Multi-scale Poisson process 
> approaches for differential expression analysis of high-throughput 
> sequencing data.* [Shim et al. (2021)][multiseq-arxiv].


## Dependency

This package depends on [the ashr R package v1.0.12](https://github.com/stephens999/ashr/releases/tag/v1.0.12). In R, you can install the ashr v1.0.12 using [devtools][devtools]:

   ```R
   install.packages("devtools")
   library(devtools)
   install_github("stephens999/ashr@v1.0.12")
   ```

## Quick Start

Follow these steps to quickly get started using multiseq.

1. In R, install the latest version of multiseq using [devtools][devtools]:

   ```R
   install.packages("devtools")
   library(devtools)
   install_github("heejungshim/multiseq")
   ```

   This will build the multiseq package *without* the vignettes. To
   build with the vignettes, do this instead:

   ```R
   install_github("heejungshim/multiseq",build_vignettes = TRUE)
   ```
   
2. Load the multiseq package:

   ```R
   library(multiseq)
   ```
   
3. To learn more, see the multiseq vignette
   (which you can also view [here][multiseq-web]):

   ```R
   vignette("multiseq")
   ```
   
[multiseq-arxiv]: https://arxiv.org/abs/2106.13634
[issues]: https://github.com/heejungshim/multiseq/issues
[gpl]: http://www.gnu.org/licenses/gpl.html
[fsf]: https://www.fsf.org
[multiseq-web]: https://heejungshim.github.io/multiseq
[devtools]: https://github.com/r-lib/devtools

