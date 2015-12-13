hgaprec: main.o hgaprec.o log.o ratings.o
	g++ -o hgaprec main.o hgaprec.o log.o ratings.o -L/usr/local/lib -L/opt/local/lib -lgsl -lgslcblas
	
main.o: main.cc env.hh hgaprec.hh log.hh
	g++ -c -std=c++11 main.cc -I. -I/usr/local/include -I/opt/local/include
	
hgaprec.o: hgaprec.cc env.hh hgaprec.hh ratings.hh gpbase.hh
	g++ -c -std=c++11 hgaprec.cc -I. -I/usr/local/include -I/opt/local/include
	
log.o: log.cc log.hh
	g++ -c -std=c++11 log.cc -I. -I/usr/local/include -I/opt/local/include
	
ratings.o: ratings.hh log.hh matrix.hh env.hh
	g++ -c -std=c++11 ratings.cc -I. -I/usr/local/include -I/opt/local/include
	
clean: 
	rm hgaprec main.o hgaprec.o log.o ratings.o
