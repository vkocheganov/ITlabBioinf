run: bayes
	./bayes
bayes: bayes.cpp
	g++ -O3 bayes.cpp  -lgsl -lgslcblas  -o bayes
