echo [DEMO] 'running weak_boundcut'
g++  weak_boundcut_weighted.cpp -o weak_boundcut_weighted.out
#for i in {1..1};do ./weak_boundcut_weighted.out ca-HepPh_SNAP.edgelist 0 0;done
for i in {1..1};do ./weak_boundcut_weighted.out ca-AstroPh_SNAP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out cit-Patents_AMINER.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out comm-EmailEnron_SNAP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out comm-WikiTalk_SNAP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out ego-twitter.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-Orkut_SNAP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-sign_epinion_snap.edgelist 2 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-sign_slashdot_snap.edgelist 2 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-SinaWeibo_NETREP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-Twitter_ASU.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-Twitter_ICWSM.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-Youtube_SNAP.edgelist 0 0;done
#for i in {1..1};do ./weak_boundcut_weighted.out AMiner-Coauthor.txt 1 0;done                             
#for i in {1..1};do ./weak_boundcut_weighted.out amazon_user_art_rate_time.edgelist 3 1;done
#for i in {1..1};do ./weak_boundcut_weighted.out nov_user_msg_time.edgelist 2 1;done
#for i in {1..1};do ./weak_boundcut_weighted.out rating-StackOverflow_KONECT.edgelist 2 1;done
#for i in {1..1};do ./weak_boundcut_weighted.out rec-YelpUserBusiness_NETREP.edgelist 1 1;done
#for i in {1..1};do ./weak_boundcut_weighted.out soc-LivejournalGrpmem2007_MPI.edgelist 0 1;done            #livejournal
#for i in {1..1};do ./weak_boundcut_weighted.out out.yahoo-song 3 1;done                                    #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.epinions-rating 3 1;done                               #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.libimseti 1 0;done                                     #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.movielens-10m_rating 3 1;done                          #konect
for i in {1..1};do ./weak_boundcut_weighted.out out.bookcrossing_rating_rating 1 1;done                    #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.librec-ciaodvd-review_ratings 1 1;done                 #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.amazon-ratings 3 1;done                               #konect
#for i in {1..1};do ./weak_boundcut_weighted.out out.wang-tripadvisor 3 1;done                              #konect
#for i in {1..1};do ./weak_boundcut_weighted.out rec-movielens-ratings.edges 3 1;done                       #networkrepository
#for i in {1..3};do ./weak_boundcut_weighted.out PP-Pathways_ppi.txt 0 0;done                               #Stanford Biomedical Network Dataset Collection
echo [DEMO] 'done!'
echo
