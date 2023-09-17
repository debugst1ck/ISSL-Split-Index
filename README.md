# Crackling

**Rapid Whole-Genome Identification of High Quality CRISPR Guide RNAs with the Crackling Method**

Jacob Bradford, Timothy Chappell, and Dimitri Perrin. The CRISPR Journal. Jun 2022.410-421. http://doi.org/10.1089/crispr.2021.0102

## Preamble

> The design of CRISPR-Cas9 guide RNAs is not trivial and is a computationally demanding task. Design tools need to identify target sequences that will maximize the likelihood of obtaining the desired cut, while minimizing off-target risk. There is a need for a tool that can meet both objectives while remaining practical to use on large genomes. 
> 
> In this study, we present Crackling, a new method that is more suitable for meeting these objectives. We test its performance on 12 genomes and on data from validation studies. Crackling maximizes guide efficiency by combining multiple scoring approaches. On experimental data, the guides it selects are better than those selected by others. It also incorporates Inverted Signature Slice Lists (ISSL) for faster off-target scoring. ISSL provides a gain of an order of magnitude in speed compared with other popular tools, such as Cas-OFFinder, Crisflash, and FlashFry, while preserving the same level of accuracy. Overall, this makes Crackling a faster and better method to design guide RNAs at scale. 
>
> Crackling is available at https://github.com/bmds-lab/Crackling under the Berkeley Software Distribution (BSD) 3-Clause license.

## Dependencies
- Any `C++` compiler, such as `g++`
