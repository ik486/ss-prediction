ss-prediction
=============

Protein Secondary Structure Prediction


This prediction technique is based on the secondary structure database provided by DSSP (http://swift.cmbi.ru.nl/gv/dssp/) and Protein Data Bank. 
From DSSP Database 56272 Unique Chains are identified for sampling and the data is provided in two files.

To predict the secondary structure run the new_prediction.py program with two arguments (no_of_sampling and the single_letter_aminoacid_sequence)
Example
 python new_prediction.py 4 CVGVFQHGKVEIIANDQGNRTTPSYVAFTDTERLIGDAAKNQVAMNPTNTVFDAKRLIGRRFDDAVVQSDMKHW
              EEECECCCCEEEHHCCCCCCCCCCEEEECCCHHHHCCHHHCHCHCCCCCCEECCCHHHCCCCCCHHHHCHHHHH
              EEEEECCCCEEEHHCCCCCCCCCCHHEECCCHHHHCCHHCCHHHCCCCCCEECHHHHCCCCCCCHHHCCCHHHH
              CCCCHCHCCEEEHHCCCCCCCCCCCEEECCCHHHHCCHHCCHHHCCCCCCEECCHHHCCCCCCCHHHHHHHHHH
              EEEEEECCEEEECCCCCCCCCEECCEEECCCCEEECHHHHCCCCCCCCCEECCCCHCCCCCCCCHHHHHHHCCC
In this case the number of sampling is four so the program samples  56272 protein sequence in to four groups. Secondary structure will be predicted using each group of data separately and program generates four secondary structures. Best sequence has to be found out using other techniques, which are not discussed here.

No of sampling decides the accuracy of the prediction. It is found that if we use a sampling of 512 almost all the 56272 aminoacid sequence produce a result of more than 95% accuracy. That is out of 512 secondary structure sequence produced accuracy of the one secondary structure sequence will be more than 95%.
