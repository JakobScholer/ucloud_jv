# Running program
Individual functions can be called through main.py with arguments, fx:
```
  python3 main.py ec
```
```
  python3 main.py img_all_stringfiles
```
```
  python3 main.py img_stringfile
```
```
  python3 main.py make_cut_dag
```
```
  python3 main.py smiles_to_reactions_bb
```
```
  python3 main.py smiles_to_reactions_nb
```
```
  python3 main.py gml
```

# Running tests:
All tests can be run with the command
```
  python3 -m unittest discover test
```
Or a single test with the command
```
  python3 -m unittest test/[filename.py]
```
