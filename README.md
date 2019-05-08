# CNN Scoring for Flexible Docking

## Abstract

> Molecular docking—the prediction of binding modes and binding affinity of a molecule to a target of known structure—is a great computational tool for structure-based drug design. However, docking scoring functions are mostly empirical or knowledge-based and the flexibility of the receptor is completely neglected in most docking studies. Recent advances in the field showed that scoring functions can be effectively learnt by convolutional neural networks (CNNs). Here we want to build on top of these findings and develop a CNN scoring function for flexible docking by extending the capabilities of `gnina`—a state-of-the-art deep learning framework for molecular docking—and by building an high-quality training dataset for flexible docking.