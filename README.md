# AttachmentTrees
This is the script I wrote in summer 2024 to help my research team better understand root finding algorithms on attachment trees.

## The attachment tree model
We work with infinite sequences of random trees (T_n)_{n >= 1} such that T_1 has vertex set V(T_1) = {0} and T_n is a tree with vertex set V(T_n) = {0,...,n-1} and edge set consisting of E(T_n) = E(T_{n-1})U{(n-1,v_{n-1})} where v_{n-1} is sampled from some distribution D_{n-1} on {0,...,n-2}. The vertex 0 is called the root of every tree T_n. Knowing the joint distribution D = (D_0,D_1,...) gives us all information about (T_n)_{n >= 1}, so we may say that (T_n)_{n >= 1} is generated through attachment process D.

We are mostly interested in the special case of Uniform Attachment Trees, where the distributions D_{n-1} are independent and uniformly distributed. In that case, we that T_n ~ UA(n).

## The problem I study
We are given a tree T of size n that has been generated through some known attachment process D, but to which we removed the vertices labels. Our goal is to find the root of T with high probability. More precisely, for any number a > 0, we want to find a subset of V(T) of minimal size such that the probability that the root is in V(T) is at least 1 - a.
