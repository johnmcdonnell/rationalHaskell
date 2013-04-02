
anderson: 
	ghc -O1 runanderson.hs

docs:
	haddock -o docs --html runanderson.hs

clean:
	rm *.o
	rm *.hi
	rm runanderson
