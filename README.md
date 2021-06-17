# BigChem_LMW

## Intro

Drug discovery and development are very expensive and time-consuming for the pharmaceutical industry. There is a need to find new computational methods that improve the efficacy and efficiency of drug findings. Virtual screening (VS) is a very effective and inexpensive in silico computational method for drug design. Recently, the accessible chemical space has grown tremendously, reaching several billion compounds. This means that, for the first time, VS can access many more compounds than even the largest experimental screening (>1000-fold more!). However, this increase in library size means that the methods used until now need to change to adapt to this new chemical space.

**Methods:** We have developed a computational pipeline to transform how VS is carried out, enabling efficient and high-quality exploration of these new and vast chemical collections. We first perform an exhaustive exploration of the fragment-like chemical space using a hierarchical strategy based on the successive pharmacophore-guided molecular docking, unsupervised clustering techniques, MM-GBSA, and Dynamic Undocking. Then, the selected fragments are grown into larger molecules using a fragment evolution platform developed in the lab. In this way, we only need to explore the regions of the chemical space that show the most potential for our target, instead of an exhaustive exploration of the whole space.

## Project dependencies

- pyMDMix: http://mdmix.sourceforge.net/
- rDock: http://rdock.sourceforge.net/
- AMBER: https://ambermd.org/
- Schrodinger Prime: https://schrodinger.com/
- MOE: https://www.chemcomp.com/Products.htm
- PyMol: https://pymol.org/

## Python and R requirements
- Python_Requirements.txt
- R_requirements.txt

