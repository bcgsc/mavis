
sphinx-apidoc -f -o source/ ./../structural_variant/ --separate
sphinx-apidoc -f -o source/ ./../ --separate

make html
