echo "do you want to update git?"
read answer
if [ $answer = "y" ]; then
read message
git add *
git commit -m $message
git push
fi
echo "do you want to update mcmp3?"
if [ $answer = "y" ]; then
scp -r * mcmp3@pascal0.hpc.unimo.it:/HOME/mcmp3/AndreaPintus/tersoff_MC/
fi