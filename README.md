# MMEDIA
Simulation files for the paper "A message passing approach for decision fusion of hidden-markov observations in the presence of synchronized attacks" published at Int. Conf. Advances in Multimedia (MMEDIA), Special Track on Models and Algorithms for Spatially and Temporally Correlated Data (STCD).

Authors: Andrea Abrardo, Mauro Barni, Kassem Kallas, Benedetta Tondi

Abstract:
We consider a setup in which a Fusion Center (FC) makes a binary decision on the sequence of system states by relying on local observations provided by both honest and byzantine nodes, ie, nodes that deliberately alter the result of the local decision to induce an error at the fusion center. In this setting, we assume a Markovian information model for the status with a given transition probability that can be perfectly estimated at the FC. Hence, we consider an attacking strategy where the byzantine nodes can coordinate their attacks by producing correlated reports, with the aim of mimicking the behavior of the original information and at the same time minimizing the information conveyed to the FC about the sequence of states. In this scenario, we derive a nearly-optimal fusion scheme based on message passing (MP) and factor graphs. Experimental results show that, although the proposed detector is able to mitigate the effect of Byzantines, the coordination of the efforts is very harmful and significantly impairs the detection performance.

Files:
Main_sim_turbo_synch_biz_improved.m: to setup the simulation parameters and run the experiments.
The rest of the files are used to calculate the messages exchanged over the factor graph and setup the Markovian observations under synchronized and unsynchronized attacks.
The synchronized attacks files can be distinguished by the extension "synch" in the file name.

The paper can be found at: https://www.academia.edu/download/56108572/MMEDIA_final_article_published.pdf

For any inquiries please contact Dr. Kassem Kallas at k_kallas@hotmail.com
