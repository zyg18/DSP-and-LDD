echo [DEMO] 'running weak_boundcut'
g++  pruning.cpp -o pruning.out
#for i in {1..1};do ./pruning.out ca-HepPh_SNAP.edgelist 0 0;done
#for i in {1..1};do ./pruning.out comm-EmailEnron_SNAP.edgelist 0 0;done
for i in {1..1};do ./pruning.out ca-AstroPh_SNAP.edgelist 0 0;done
#for i in {1..1};do ./pruning.out PP-Pathways_ppi.txt 0 0;done                               #Stanford Biomedical Network Dataset Collection
#for i in {1..1};do ./pruning.out soc-Twitter_ICWSM.edgelist 0 0;done
#for i in {1..1};do ./pruning.out soc-sign_slashdot_snap.edgelist 2 0;done
#for i in {1..1};do ./pruning.out rating-StackOverflow_KONECT.edgelist 2 1;done
#for i in {1..1};do ./pruning.out soc-sign_epinion_snap.edgelist 2 0;done
#for i in {1..1};do ./pruning.out ego-twitter.edgelist 0 0;done
#for i in {1..1};do ./pruning.out soc-Youtube_SNAP.edgelist 0 0;done
#for i in {1..1};do ./pruning.out comm-WikiTalk_SNAP.edgelist 0 0;done
#for i in {1..1};do ./pruning.out nov_user_msg_time.edgelist 2 1;done
#for i in {1..1};do ./pruning.out cit-Patents_AMINER.edgelist 0 0;done
#for i in {1..1};do ./pruning.out soc-Twitter_ASU.edgelist 0 0;done
#for i in {1..1};do ./pruning.out soc-LivejournalGrpmem2007_MPI.edgelist 0 1;done            #livejournal
#for i in {1..1};do ./pruning.out soc-Orkut_SNAP.edgelist 0 0;done
#for i in {1..1};do ./pruning.out soc-SinaWeibo_NETREP.edgelist 0 0;done

#for i in {1..1};do ./pruning.out out.wang-tripadvisor 3 1;done                              #konect
#for i in {1..1};do ./pruning.out rec-YelpUserBusiness_NETREP.edgelist 1 1;done
for i in {1..1};do ./pruning.out out.bookcrossing_rating_rating 1 1;done                    #konect
#for i in {1..1};do ./pruning.out out.librec-ciaodvd-review_ratings 1 1;done                 #konect
#for i in {1..1};do ./pruning.out out.movielens-10m_rating 3 1;done                          #konect
#for i in {1..1};do ./pruning.out out.epinions-rating 3 1;done                               #konect
#for i in {1..1};do ./pruning.out out.libimseti 1 0;done                                     #konect
#for i in {1..1};do ./pruning.out rec-movielens-ratings.edges 3 1;done                       #networkrepository
#for i in {1..1};do ./pruning.out out.yahoo-song 3 1;done                                    #konect
echo [DEMO] 'done!'
echo
