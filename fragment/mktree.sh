root -l -b -q 'mktree_evt.C(13,13,0,"24Mar")' > 13empty.log &
root -l -b -q 'mktree_evt.C(13,13,1,"24Mar")' > 13carbon.log &
root -l -b -q 'mktree_evt.C(13,13,2,"24Mar")' > 13ch2.log &
root -l -b -q 'mktree_evt.C(11,122,0,"24Mar")' > 12empty.log &
root -l -b -q 'mktree_evt.C(11,122,1,"24Mar")' > 12carbon.log &
root -l -b -q 'mktree_evt.C(11,122,2,"24Mar")' > 12ch2.log &
watch 'tail 1*log'
