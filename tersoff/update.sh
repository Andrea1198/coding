echo "do you want to update git? (y/n)"
read answer
if [ $answer = "y" ]; then
echo "insert a message for git"
read message
git add *
git commit -m $message
git push
fi
echo "do you want to update mcmp3? (y/n)"
read answer
# if [ $answer = "y" || $answer = "Y" || $answer = "yes" || $answer = "YES"]; then
if [ $answer = "y" ]; then
scp -r * mcmp3@pascal0.hpc.unimo.it:/HOME/mcmp3/AndreaPintus/tersoff_MC/
fi
