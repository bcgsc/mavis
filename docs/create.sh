
sphinx-apidoc -f -P -o source/ ./../structural_variant/  --separate
sphinx-apidoc -f -P -o source/ ./../tools --separate


cp ./../README.rst source/
system=$(uname -s)

for x in source/structural_variant.*rst
do
    echo $x;
    if [ "$system" == "Darwin" ]; then # mac osx
        sed -i '.bk' 's/:show-inheritance:/:special-members: __and__, __or__, __xor__, __len__, __sub__, __add__/g' $x;
        rm "$x.bk"
    else # for others assume linux
        sed -i 's/:show-inheritance:/:special-members: __and__, __or__, __xor__, __len__, __sub__, __add__/g' $x; # LINUX
    fi
done
make html
