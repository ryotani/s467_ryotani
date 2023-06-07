root -l -b -q 'mktree_evt.C(13,13,0,"23May")' > 13empty.log &
root -l -b -q 'mktree_evt.C(13,13,1,"23May")' > 13carbon.log &
root -l -b -q 'mktree_evt.C(13,13,2,"23May")' > 13ch2.log &
root -l -b -q 'mktree_evt.C(11,122,0,"23May")' > 12empty.log &
root -l -b -q 'mktree_evt.C(11,122,1,"23May")' > 12carbon.log &
root -l -b -q 'mktree_evt.C(11,122,2,"23May")' > 12ch2.log &
watch 'tail 1*log'
