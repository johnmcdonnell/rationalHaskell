
bias:
	ghc -O1 randombias.hs

anderson: 
	ghc -O1 testanderson.hs

clean:
	rm *.o
	rm *.hi
	rm testanderson
