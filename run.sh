for input in input/g*.txt; do
    echo ""
    echo "Copy an instance to maxcut.in file."
    echo "cp input/g1.txt maxcut.in"
    cp $input maxcut.in

    echo ""
    echo "make all"
    make all

    echo ""
    echo "make run"
    make run

    echo ""
    echo " * 'make clean' before submit or re-run"
    echo " * Remove irrelevant print functions before submit."
    echo " * Your program must end before time limit."
    echo ""
done