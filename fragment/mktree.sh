root -l -b -q 'mktree_evt.C(13,13,0,"20Jul")' > 13empty.log &
root -l -b -q 'mktree_evt.C(13,13,1,"20Jul")' > 13carbon.log &
root -l -b -q 'mktree_evt.C(13,13,2,"20Jul")' > 13ch2.log &
root -l -b -q 'mktree_evt.C(11,122,0,"20Jul")' > 12empty.log &
root -l -b -q 'mktree_evt.C(11,122,1,"20Jul")' > 12carbon.log &
root -l -b -q 'mktree_evt.C(11,122,2,"20Jul")' > 12ch2.log &
watch 'tail 1*log'
