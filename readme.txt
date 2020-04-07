%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This is  a matlab package for  robust decoding from one bit measurements via LQ-constrained Least Squares. 
The main reference is the package about the one-bit compressive sensing via L1-regularized Least Squares. 
This package is  maintained by Cui Jia
Email: cjia80@whu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rec_selec_vote.m to verify the efficiency of the majority vote parameter selection rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Other files reproduces all the numerical results in the paper 
 Fig.1 by demo_1bit_LQhigh .m 
 Fig.2 by rec_prob_s.m, rec_prob_sigma.m and  rec_prob_q.m
 Fig.3 by rec_s_compare.m
Fig.4 by rec_sigma_compare.m and  rec_sign_compare.m
Fig.5 by rec_sigma_compare.m  and  rec_sign_compare.m
Tab.1 by Table1_compare.m, 
 Fig.6 and Tab.2 by onebit_1d.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
If you have any questions or find any buggs please contact 
cjia80@whu.edu.cn and  ly13781573969@163.com. Thank you.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Methods are:
1. wpdasc: download from http://faculty.zuel.edu.cn/tjyjxxy/jyl/list.htm.
      which is an implementation by Cui Jia according to this paper.
2. LinProj: download from http://faculty.zuel.edu.cn/tjyjxxy/jyl/list.htm/
   which is an implementation by Yuling Jiao according to the paper of 
   Roman Vershynin  "Estimation in high dimensions: a geometric perspective,
   in Sampling theory, a renaissance, Springer, 2015, pp. 3¨C66."
    
3. BIHT_1: download from http://perso.uclouvain.be/laurent.jacques/index.
   php/Main/BIHTDemo, L. Jacques, J.N. Laska, P.T. Boufounos, and 
   R.G. Baraniuk, ¡°Robust 1-bit compressive sensing via binary stable 
   embeddings of sparse vectors,¡± IEEE Transactions on Information Theory,
    vol. 59, no. 4, pp. 2082--2102, 2013.

4. BIHT_AOP_flip: download from http://www.esat.kuleuven.be/stadius/ADB/
   huang/downloads/1bitCSLab.zip, M. Yan, Y. Yang, and S. Osher, "Robust 
   1-bit compressive sensing using adaptive outlier pursuit," IEEE 
   Transactions on Signal Processing, vol. 60, no. 7, pp. 3868--3875, 2012.

5. PIHT_AOP_flip: download from http://www.esat.kuleuven.be/stadius/ADB/
   huang/downloads/1bitCSLab.zip, X. Huang, L. Shi, M. Yan, J.A.K. Suykens,
   Pinball loss minimization for one-bit compressive sensing. Internal 
   Report 15-76, ESAT-SISTA, KU Leuven. 

6.pdasc: download from http://faculty.zuel.edu.cn/tjyjxxy/jyl/list.htm.

7. mcp_1bit£ºreference: "Xiaolin Huang and Ming Yan. Nonconvex penalties with analytical solutions for one-bit compressive
   sensing. Signal Processing, 144:341¨C351, 2018".
   which is an implementation by Cui Jia taken the sparsity prior knowledge into consideration.

8. passive_1bit: download http://www.esat.kuleuven.be/stadius/ADB/huang/downloads/1bitCSLab.zip.
   and implementation according  "Xiaolin Huang and Ming Yan. Nonconvex penalties with analytical solutions for one-bit compressive
   sensing. Signal Processing, 144:341¨C351, 2018".
   which is an implementation by Cui Jia taken the sparsity prior knowledge into consideration.