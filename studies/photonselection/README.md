# Photon Selection

The purpose of this study is to determine the optimal means to choose candidate neutrino vertices for the photon selection.

We aim to develop a selection module that can then be run as part of the main gen2ntuple program.

To start, we will bias the selected neutrino candidate towards events

* with a shower prong
* is consistent with the in-time flash
* has the largest total ionization

Our study operates on the output of the lantern reco, and therefore the contents of the KPSRecoManagerTree.
We augment this information with the flash prediction.

Our first study is to compile selection variables that might help us discriminate true versus false photon vertices.

So one key task is to define a truth selection. Our target events:

* has a primary photon that deposits enough energy in the detector
* we can define a primary photon as one that originates near the vertex. they are almost always the result of a neutral pion decay.
* we can also indicate what other particles from the neutrino interaction leave behind ionization in the TPC

