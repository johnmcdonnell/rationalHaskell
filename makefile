
anderson: 
	ghc -O1 testanderson.hs

docs:
	haddock -o docs --html testanderson.hs

clean:
	rm *.o
	rm *.hi
	rm testanderson
