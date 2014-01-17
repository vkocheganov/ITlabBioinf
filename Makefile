run: NaiveBayes
	./NaiveBayes
NaiveBayes: NaiveBayes.cpp
	g++ -O3 NaiveBayes.cpp  -lgsl -lgslcblas  -o NaiveBayes
