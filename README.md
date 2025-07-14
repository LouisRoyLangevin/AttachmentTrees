# AttachmentTrees
This is the script I wrote in summer 2024 to help my research team better understand root finding algorithms on attachment trees.

## The attachment tree model
We work with infinite sequences of random trees (T_n)_{n >= 1} such that T_1 has vertex set V(T_1) = {0} and T_n is a tree with vertex set V(T_n) = {0,...,n-1} and edge set consisting of E(T_n) = E(T_{n-1})U{(n-1,v_{n-1})} where v_{n-1} is sampled from some distribution D_{n-1} on {0,...,n-2}. The vertex 0 is called the root of every tree T_n. Knowing the joint distribution D = (D_0,D_1,...) gives us all information about (T_n)_{n >= 1}, so we may say that (T_n)_{n >= 1} is generated through attachment process D.

We are mostly interested in the special case of Uniform Attachment trees, where the distributions D_{n-1} are independent and uniformly distributed. In that case, we that T_n ~ UA(n).

## The problem I study
We are given a tree T of size n that has been generated through some known attachment process D, but to which we removed the vertices labels. Our goal is to find the root of T with high probability. More precisely, for any number a > 0, we want to find a subset of V(T) of minimal size N = N(T,a) such that the probability that the root is in V(T) is at least 1 - a.

From now on, we are working with Uniform Attachment trees. It turns out that a paper published by Bubeck,Devroye,Lugosi in 2015 showed that for any number a > 0, there exists an (minimal) integer K = K(a) such that N(T,a) <= K for any such tree T. In fact, this same paper showed that there exists a universal constant c > 0 such that K >= exp(c + c*\sqrt{log(1/a)}). 

## What this script is about

In 2024, my research team and I proved that there also exists a universal constant c' > 0 such that K <= exp(c' + c'*\sqrt{log(1/a)}), tightening the bound up to a constant factor. This solved one of the most important problems in the subfield of root-finding algorithms in the last 10 years. However, it still remains to find the minimal possible value of c'.

This is a very hard task, so to give ourselves an idea, I decided to write this script which binary searches through the possible values of c' and approximates the minimal value. I wrote this script under the assumption that K(a) "reaches its asymptotic behaviour" before a = 2e-4 (because putting a > 0 too small would make the program very slow).

## For more information

I refer to our paper ["OPTIMAL ROOT RECOVERY FOR UNIFORM ATTACHMENT TREES AND d-REGULAR GROWING TREES"](https://arxiv.org/abs/2411.18614) published in 2024 on Arxiv. This paper contains more information about the process used to find the vertices with highest probability of being the root. This includes the definition of the logPhi function that you can find in the script.
