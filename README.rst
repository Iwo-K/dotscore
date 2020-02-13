.. role:: pyth(code)
  :language: python

DoTscore
========

dotscore is a python module created to enable easy computation of DoT-scores (Direction of Transition scores) as presented in Kucinski et al. 2020 (link).
DoT-score is method aiding interpretation of transcriptional changes (e.g. lists of differentially expressed genes) using scRNA-Seq landscapes as a reference. The module is built on top of the scanpy module (https://scanpy.readthedocs.io/en/stable/)

DoT-score concept
-----------------

In a scRNA-Seq landscape gene expression values can be scaled to a chosen point of origin. This creates a set of vectors connecting each cell to the point of origin - cell state vectors in the multidimensional gene expression space. 
Now consider a treatment (gene knockout, chemical perturbation etc.) to a given cell population followed by a transcriptomic readout. The changes in gene expression (e.g. log2(Fold Change) from differential expression) are also a vector - the positive and negative values indicating shifts in the many dimensions of gene expression space.

DoT-score computes how well the two vectors (cell state change and treatment change) are aligned and enable visualisation on the scRNA-Seq landscape.
A positive DoT-score for a cell on the landscape indicates that the treatment causes cells to shift towards that state (starting from the point of origin), while a negative DoT-score indicates a shift away from that state. We provide a diagram below to illustrate this concept and a concrete example here: **FILL IN....**


.. figure:: images/DoT_diagram.png
   :height: 400px
   :scale: 25 %
   :align: center


More specifically, DoT score (**s**) calculates the dot product (proportional to the angle) between the treatment vector and the vector of gene expression (scaled) for each cell on the landscape: **s** = **X** **v**, where **X** is a matrix cells x genes with scaled expression values and **v** is a vector of weights (e.g. log2(Fold Change) for each gene). 
To provide a measure of statistical significance, we simulate the DoT-score by randomly choosing weights and genes and computing a z-score.

Installation
------------

Python > 3.4 and pip are required. To install the package:

1. Clone the repository:

.. code-block:: text

    git clone https://github.com/Iwo-K/dotscore

2. Install the dependencies

.. code-block:: text

    pip install -r ./dotscore/requirements.txt

3. Install the package

.. code-block:: text

    pip install -e ./dotscore/

Usage
-----

Dotscore is based on the AnnData/scanpy framework, thus AnnData objects are starting points for most functions. Module exports the following functions:

- custom_scale
- get_DoTscore
- 

To compute the DoT-score for cell of a given 


make_app() accepts the following arguments:
  - :pyth:`adata` - an AnnData object

