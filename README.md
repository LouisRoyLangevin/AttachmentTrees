# AttachmentTrees
This is the script I wrote in summer 2024 to help my research team better understand root finding algorithms on attachment trees.

## The attachment tree model
We work with infinite sequences of random trees (T_n) such that T_1 has vertex set V(T_1) = {0} and T_n is a tree with vertex set V(T_n) = {0,...,n-1} and edge set consisting of E(T_n) = E(T_{n-1})U{(n-1,v_{n-1})} where v_{n-1} is sampled from some distribution D_{n-1} on {0,...,n-2}.

We are mostly interested in the special case of Uniform Attachment Trees, where the distributions D_{n-1} are independent and uniformly distributed. In that case, we that T_n ~ UA(n).
