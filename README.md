# Grid Switch Stabilisation

This repository contains the implementation of a numerical procedure meant to compute a switch stabilizing control law for a power distribution network.

## Non-distributed controller

Initially, a stabilizing baseline controller is realised in a switch stabilizing manner but with no sparsity guarantees. Details can be found in the files inside the folder *5_nodes_network*.

## Distributed controller

The files inside *5_nodes_distributed* contain the numerical procedures utilised to compute a switch stabilizing controller with a pre-imposed sparsity pattern. This is done with the help of a sequential program capable of solving BLMIs.

## Distributed controller - IMC

Finally, the above-mentioned controller is augmented using the internal model control principle with the aim of ensuring rejection of step-like disturbance signals.

---

More details on the theoretical aspects can be found in the linked article: LINK