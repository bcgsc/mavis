
sphinx-apidoc -f -P -o source/ ./../structural_variant/ --separate
sphinx-apidoc -f -P -o source/ ./../ --separate
cp ./../README.md source/
for x in source/structural_variant.*rst
do
    echo $x;
    sed -i '.bk' 's/:show-inheritance:/:special-members: __and__, __or__, __xor__, __len__, __sub__, __add__/g' $x;
done
make html
