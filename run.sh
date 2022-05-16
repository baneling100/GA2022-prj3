echo ""
echo "Copy an instance to maxcut.in file."
echo "cp input/g1.txt maxcut.in"
cp input/g1.txt maxcut.in

echo ""
echo "make all"
make all

echo ""
echo "make run"
timeout 180 make run 2>&1

echo ""
echo " * 'make clean' before submit or re-run"
echo " * Remove irrelevant print functions before submit."
echo " * Your program must end before time limit."
echo ""
